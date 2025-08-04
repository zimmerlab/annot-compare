package com.github.zimmerlab.gtfcompare.utils;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.model.TranscriptPair;
import com.github.zimmerlab.gtfcompare.model.comparison.TranscriptComparisonResult;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

enum FeatureType {UTR5, CDS, UTR3, INTRON, START_CODON, STOP_CODON}

public class OverlappingTranscripts {
    public static List<TranscriptPair> map(GtfFile targetGtfFile, GtfFile queryGtfFile) {

        // interval trees
        Map<String, IntervalTreeMap<List<TranscriptFeature>>> targetTrees = buildIntervalTrees(targetGtfFile);
        Map<String, IntervalTreeMap<List<TranscriptFeature>>> queryTrees = buildIntervalTrees(queryGtfFile);

        // find overlaps
        List<TranscriptPair> allPairs = findOverlaps(targetTrees, queryTrees).toList();

        // cluster loci
        Map<TranscriptFeature, List<TranscriptPair>> clusters = clusterOverlappingPairs(allPairs);

        // best hit per loci
        var finalMapping = new ConcurrentHashMap<TranscriptFeature, TranscriptFeature>();
        clusters.values().parallelStream().forEach(locusPairs -> {
            Set<TranscriptFeature> targets = locusPairs.stream().map(TranscriptPair::getTargetTranscript).collect(Collectors.toCollection(LinkedHashSet::new));
            Set<TranscriptFeature> queries = locusPairs.stream().map(TranscriptPair::getQueryTranscript).collect(Collectors.toCollection(LinkedHashSet::new));

            List<TranscriptFeature> clique = Stream.concat(targets.stream(), queries.stream()).collect(Collectors.toList());
            var vectors = buildModelVectorsForClique(clique);

            Map<TranscriptFeature, TranscriptFeature> mapping = mapReciprocalBestHits(targets, queries, vectors);
            finalMapping.putAll(mapping);
        });

        //debug(allPairs, targetGtfFile, queryGtfFile, finalMapping);

        return finalMapping.entrySet().stream().map(entry -> new TranscriptPair(entry.getKey(), entry.getValue(), new TranscriptComparisonResult())).toList();
    }

    private static void debug(List<TranscriptPair> allPairs, GtfFile targetGtfFile, GtfFile queryGtfFile, Map<TranscriptFeature, TranscriptFeature> finalMapping) {
        Map<TranscriptFeature, List<TranscriptPair>> overlapsByTarget =
                allPairs.stream()
                        .collect(Collectors.groupingBy(
                                TranscriptPair::getTargetTranscript,
                                Collectors.toList()
                        ));

        Map<TranscriptFeature, List<TranscriptPair>> overlapsByQuery =
                allPairs.stream()
                        .collect(Collectors.groupingBy(
                                TranscriptPair::getQueryTranscript,
                                Collectors.toList()
                        ));
        Set<TranscriptFeature> allTargets = getAllTranscripts(targetGtfFile);
        Set<TranscriptFeature> allQueries = getAllTranscripts(queryGtfFile);

        Set<TranscriptFeature> unmappedTargets = new LinkedHashSet<>(allTargets);
        unmappedTargets.removeAll(finalMapping.keySet());

        Set<TranscriptFeature> noOverlapTargets = unmappedTargets.stream()
                .filter(t -> !overlapsByTarget.containsKey(t))
                .collect(Collectors.toCollection(LinkedHashSet::new));

        Set<TranscriptFeature> noRecipBestHitTargets = unmappedTargets.stream()
                .filter(t -> overlapsByTarget.containsKey(t))
                .collect(Collectors.toCollection(LinkedHashSet::new));

        Set<TranscriptFeature> unmappedQueries = new LinkedHashSet<>(allQueries);
        unmappedQueries.removeAll(new HashSet<>(finalMapping.values()));

        Set<TranscriptFeature> noOverlapQueries = unmappedQueries.stream()
                .filter(q -> !overlapsByQuery.containsKey(q))
                .collect(Collectors.toCollection(LinkedHashSet::new));

        Set<TranscriptFeature> noRecipBestHitQueries = unmappedQueries.stream()
                .filter(q -> overlapsByQuery.containsKey(q))
                .collect(Collectors.toCollection(LinkedHashSet::new));

        System.out.println("=== Targets ohne jegliches Overlap ===");
        noOverlapTargets.forEach(System.out::println);

        System.out.println("=== Targets mit Overlap, aber ohne reziproken Best-Hit ===");
        noRecipBestHitTargets.forEach(System.out::println);

        System.out.println("=== Queries ohne jegliches Overlap ===");
        noOverlapQueries.forEach(System.out::println);

        System.out.println("=== Queries mit Overlap, aber ohne reziproken Best-Hit ===");
        noRecipBestHitQueries.forEach(System.out::println);
    }

    private static Set<TranscriptFeature> getAllTranscripts(GtfFile gtfFile) {
        return gtfFile.getAllGeneFeatureIds().stream()
                .flatMap(geneId ->
                        gtfFile.getGeneFeature(geneId)
                                .getTranscripts().stream()
                )
                .collect(Collectors.toCollection(LinkedHashSet::new));
    }


    private static Map<String, IntervalTreeMap<List<TranscriptFeature>>> buildIntervalTrees(GtfFile gtfFile) {
        var trees = new HashMap<String, IntervalTreeMap<List<TranscriptFeature>>>();
        Function<TranscriptFeature, String> keyOf = t -> t.getBaseData().getContig();
        extractTranscripts(gtfFile, trees, keyOf);
        return trees;
    }

    private static Map<TranscriptFeature, List<TranscriptPair>> clusterOverlappingPairs(List<TranscriptPair> allPairs) {
        OverlappingTranscripts.DisjointSet<TranscriptFeature> ds = new OverlappingTranscripts.DisjointSet<>();
        allPairs.forEach(p -> ds.union(p.getTargetTranscript(), p.getQueryTranscript()));

        return allPairs.stream().collect(Collectors.groupingBy(p -> ds.find(p.getTargetTranscript()), LinkedHashMap::new, Collectors.toList()));
    }

    private static Map<TranscriptFeature, TranscriptFeature> mapReciprocalBestHits(Set<TranscriptFeature> targets, Set<TranscriptFeature> queries, Map<TranscriptFeature, EnumMap<FeatureType, BitSet>> vectors) {
        var bestQPerT = new HashMap<TranscriptFeature, TranscriptFeature>();
        var bestScoreQ = new HashMap<TranscriptFeature, Double>();
        for (TranscriptFeature t : targets) {
            bestScoreQ.put(t, -1.0);
            var vecT = vectors.get(t);
            for (TranscriptFeature q : queries) {
                double score = jaccardSimilarity(vecT, vectors.get(q));
                if (score > bestScoreQ.get(t)) {
                    bestScoreQ.put(t, score);
                    bestQPerT.put(t, q);
                    if (score >= 1.0) break;
                }
            }
        }

        var bestTPerQ = new HashMap<TranscriptFeature, TranscriptFeature>();
        var bestScoreT = new HashMap<TranscriptFeature, Double>();
        for (TranscriptFeature q : queries) {
            bestScoreT.put(q, -1.0);
            var vecQ = vectors.get(q);
            for (TranscriptFeature t : targets) {
                double score = jaccardSimilarity(vectors.get(t), vecQ);
                if (score > bestScoreT.get(q)) {
                    bestScoreT.put(q, score);
                    bestTPerQ.put(q, t);
                }
            }
        }

        var mapping = new LinkedHashMap<TranscriptFeature, TranscriptFeature>();
        for (var e : bestQPerT.entrySet()) {
            TranscriptFeature t = e.getKey();
            TranscriptFeature q = e.getValue();
            if (t.equals(bestTPerQ.get(q))) {
                mapping.put(t, q);
            }
        }
        return mapping;
    }

    private static Stream<TranscriptPair> findOverlaps(Map<String, IntervalTreeMap<List<TranscriptFeature>>> targetTrees, Map<String, IntervalTreeMap<List<TranscriptFeature>>> queryTrees) {

        return targetTrees.entrySet().stream().flatMap(entry -> {
            String contig = entry.getKey();
            var tTree = entry.getValue();
            var qTree = queryTrees.get(contig);
            if (qTree == null) return Stream.empty();

            return tTree.entrySet().stream().flatMap(e -> {
                var iv = e.getKey();
                var tList = e.getValue();
                return qTree.getOverlapping(iv).stream().flatMap(List::stream).flatMap(q -> tList.stream().map(t -> new TranscriptPair(t, q)));
            });
        });
    }


    private static double jaccardSimilarity(EnumMap<FeatureType, BitSet> a, EnumMap<FeatureType, BitSet> b) {
        int inter = 0, uni = 0;
        for (FeatureType ft : FeatureType.values()) {
            BitSet ai = (BitSet) a.get(ft).clone();
            ai.and(b.get(ft));
            BitSet au = (BitSet) a.get(ft).clone();
            au.or(b.get(ft));
            inter += ai.cardinality();
            uni += au.cardinality();
        }
        return uni > 0 ? (double) inter / uni : 0.0;
    }


    private static void extractTranscripts(GtfFile gtfFile, Map<String, IntervalTreeMap<List<TranscriptFeature>>> trees, Function<TranscriptFeature, String> keyOf) {
        for (var geneId : gtfFile.getAllGeneFeatureIds()) {
            for (var tf : gtfFile.getGeneFeature(geneId).getTranscripts()) {
                String key = keyOf.apply(tf);
                var bd = tf.getBaseData();
                var iv = new Interval(bd.getContig(), bd.getStart(), bd.getEnd());

                trees.computeIfAbsent(key, k -> new IntervalTreeMap<>()).compute(iv, (interval, list) -> {
                    if (list == null) list = new ArrayList<>();
                    list.add(tf);
                    return list;
                });
            }
        }
    }


    private static class DisjointSet<T> {
        private final Map<T, T> parent = new HashMap<>();

        T find(T x) {
            parent.putIfAbsent(x, x);
            if (!parent.get(x).equals(x)) {
                parent.put(x, find(parent.get(x)));
            }
            return parent.get(x);
        }

        void union(T a, T b) {
            T ra = find(a), rb = find(b);
            if (!ra.equals(rb)) parent.put(ra, rb);
        }
    }

    private static Map<TranscriptFeature, EnumMap<FeatureType, BitSet>> buildModelVectorsForClique(List<TranscriptFeature> clique) {
        int locusStart = clique.stream().mapToInt(t -> t.getBaseData().getStart()).min().orElseThrow();
        int locusEnd = clique.stream().mapToInt(t -> t.getBaseData().getEnd()).max().orElseThrow();
        int length = locusEnd - locusStart + 1;

        Map<TranscriptFeature, EnumMap<FeatureType, BitSet>> all = new LinkedHashMap<>();

        for (var tf : clique) {
            EnumMap<FeatureType, BitSet> map = new EnumMap<>(FeatureType.class);
            for (FeatureType ft : FeatureType.values()) {
                map.put(ft, new BitSet(length));
            }

            for (var seg : tf.getFeatures()) {
                FeatureType ft = getFeatureType(seg);

                if (ft == null) continue;

                int s = seg.getBaseData().getStart() - locusStart;
                int e = seg.getBaseData().getEnd() - locusStart;
                map.get(ft).set(s, e + 1);
            }

            all.put(tf, map);
        }

        return all;
    }

    private static FeatureType getFeatureType(GtfFeature seg) {
        String type = seg.getBaseData().getType();
        FeatureType ft = null;

        if (Constants.UTR5.equals(type)) {
            ft = FeatureType.UTR5;
        } else if (Constants.CDS.equals(type)) {
            ft = FeatureType.CDS;
        } else if (Constants.UTR3.equals(type)) {
            ft = FeatureType.UTR3;
        } else if (Constants.INTRON.equals(type)) {
            ft = FeatureType.INTRON;
        } else if (Constants.START_CODON.equals(type)) {
            ft = FeatureType.START_CODON;
        } else if (Constants.STOP_CODON.equals(type)) {
            ft = FeatureType.STOP_CODON;
        }
        return ft;
    }
}
