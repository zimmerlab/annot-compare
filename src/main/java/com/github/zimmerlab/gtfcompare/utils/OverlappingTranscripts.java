package com.github.zimmerlab.gtfcompare.utils;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.model.TranscriptPair;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class OverlappingTranscripts {
    public static Map<TranscriptFeature, TranscriptFeature> map(GtfFile targetGtfFile, GtfFile queryGtfFile) {

        var targetTrees = new HashMap<String, IntervalTreeMap<TranscriptFeature>>();
        var queryTrees = new HashMap<String, IntervalTreeMap<TranscriptFeature>>();
        Function<TranscriptFeature, String> keyOf = t -> t.getBaseData().getContig();
        extractTranscripts(targetGtfFile, targetTrees, keyOf);
        extractTranscripts(queryGtfFile, queryTrees, keyOf);

        var allPairs = findOverlaps(targetTrees, queryTrees).toList();

        var ds = new OverlappingTranscripts.DisjointSet<TranscriptFeature>();
        allPairs.forEach(p -> ds.union(p.getTargetTranscript(), p.getQueryTranscript()));

        var rootToPairs = allPairs.stream().collect(Collectors.groupingBy(p -> ds.find(p.getTargetTranscript()), LinkedHashMap::new, Collectors.toList()));

        var finalMapping = new LinkedHashMap<TranscriptFeature, TranscriptFeature>();
        for (List<TranscriptPair> locusPairs : rootToPairs.values()) {
            var targets = locusPairs.stream().map(TranscriptPair::getTargetTranscript).collect(Collectors.toCollection(LinkedHashSet::new));
            var queries = locusPairs.stream().map(TranscriptPair::getQueryTranscript).collect(Collectors.toCollection(LinkedHashSet::new));

            var clique = Stream.concat(targets.stream(), queries.stream()).collect(Collectors.toList());
            var vectors = buildModelVectorsForClique(clique);

            for (TranscriptFeature t : targets) {
                var vt = vectors.get(t);
                TranscriptFeature bestQ = null;
                double bestScore = -1.0;
                for (TranscriptFeature q : queries) {
                    double score = jaccardSimilarity(vt, vectors.get(q));
                    if (score > bestScore) {
                        bestScore = score;
                        bestQ = q;
                    }
                }
                if (bestQ != null) {
                    finalMapping.put(t, bestQ);
                }
            }
        }

        List<TranscriptFeature> allTargets = targetGtfFile.getAllGeneFeatureIds().stream().flatMap(geneId -> targetGtfFile.getGeneFeature(geneId).getTranscripts().stream()).toList();

        int totalTargets = allTargets.size();

        Set<TranscriptFeature> overlappedTargets = allPairs.stream().map(TranscriptPair::getTargetTranscript).collect(Collectors.toSet());
        int overlapCount = overlappedTargets.size();

        int mappedCount = finalMapping.size();
        int failedWithinOverlap = overlapCount - mappedCount;
        int noOverlap = totalTargets - overlapCount;

        System.out.printf("totalTargets = %d%n", totalTargets);
        System.out.printf("overlapCount = %d%n", overlapCount);
        System.out.printf("mappedCount = %d%n", mappedCount);
        System.out.printf("failedWithinOverlap = %d%n", failedWithinOverlap);
        System.out.printf("noOverlap = %d%n", noOverlap);


        List<TranscriptFeature> nonOverlappingTargets = allTargets.stream().filter(t -> !overlappedTargets.contains(t)).toList();

        return finalMapping;
    }

    private static Stream<TranscriptPair> findOverlaps(Map<String, IntervalTreeMap<TranscriptFeature>> targetTrees, Map<String, IntervalTreeMap<TranscriptFeature>> queryTrees) {

        return targetTrees.entrySet().stream().flatMap(entry -> {
            String contig = entry.getKey();
            IntervalTreeMap<TranscriptFeature> tTree = entry.getValue();
            IntervalTreeMap<TranscriptFeature> qTree = queryTrees.get(contig);
            if (qTree == null) return Stream.empty();
            return tTree.values().stream().flatMap(t -> {
                var td = t.getBaseData();
                Interval iv = new Interval(td.getContig(), td.getStart(), td.getEnd());
                return qTree.getOverlapping(iv).stream().map(q -> new TranscriptPair(t, q));
            });
        });
    }


    private static double jaccardSimilarity(BitSet a, BitSet b) {
        var intersection = (BitSet) a.clone();
        intersection.and(b);
        var union = (BitSet) a.clone();
        union.or(b);
        int interSize = intersection.cardinality();
        int unionSize = union.cardinality();
        return unionSize > 0 ? (double) interSize / unionSize : 0.0;
    }


    private static void extractTranscripts(GtfFile gtfFile, HashMap<String, IntervalTreeMap<TranscriptFeature>> trees, Function<TranscriptFeature, String> keyOf) {
        for (var id : gtfFile.getAllGeneFeatureIds()) {
            var g = gtfFile.getGeneFeature(id);
            var t = g.getTranscripts();

            for (var transcript : t) {
                String key = keyOf.apply(transcript);
                var baseData = transcript.getBaseData();
                trees.computeIfAbsent(key, k -> new IntervalTreeMap<>()).put(new Interval(baseData.getContig(), baseData.getStart(), baseData.getEnd()), transcript);
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

    private static Map<TranscriptFeature, BitSet> buildModelVectorsForClique(List<TranscriptFeature> clique) {
        int locusStart = clique.stream().mapToInt(t -> t.getBaseData().getStart()).min().orElseThrow();
        int locusEnd = clique.stream().mapToInt(t -> t.getBaseData().getEnd()).max().orElseThrow();
        int length = locusEnd - locusStart + 1;

        char[] base = new char[length];
        Arrays.fill(base, 'G');

        Map<TranscriptFeature, BitSet> vectors = new LinkedHashMap<>();

        for (TranscriptFeature tf : clique) {
            char[] model = Arrays.copyOf(base, length);

            var exons = tf.getFeatures("exon").stream().map(e -> new Interval(e.getBaseData().getContig(), e.getBaseData().getStart(), e.getBaseData().getEnd())).sorted(Comparator.comparingInt(Interval::getStart)).toList();
            for (int i = 0; i + 1 < exons.size(); i++) {
                int intronStart = exons.get(i).getEnd() + 1;
                int intronEnd = exons.get(i + 1).getStart() - 1;
                for (int pos = intronStart; pos <= intronEnd; pos++) {
                    model[pos - locusStart] = 'I';
                }
            }

            for (var seg : tf.getFeatures()) {
                int segStart = seg.getBaseData().getStart();
                int segEnd = seg.getBaseData().getEnd();
                String type = seg.getBaseData().getType();
                char code;
                if (Constants.UTR5.equals(type)) code = 'F';
                else if (Constants.CDS.equals(type)) code = 'C';
                else if (Constants.UTR3.equals(type)) code = 'T';
                else if (Constants.INTRON.equals(type)) code = 'I';
                else if (Constants.START_CODON.equals(type)) code = 'S';
                else if (Constants.STOP_CODON.equals(type)) code = 'E';
                else code = 'G';
                for (int pos = segStart; pos <= segEnd; pos++) {
                    model[pos - locusStart] = code;
                }
            }

            var v = new BitSet(length);
            for (int i = 0; i < length; i++) {
                if (model[i] != 'G') {
                    v.set(i);
                }
            }

            vectors.put(tf, v);
        }

        return vectors;
    }
}
