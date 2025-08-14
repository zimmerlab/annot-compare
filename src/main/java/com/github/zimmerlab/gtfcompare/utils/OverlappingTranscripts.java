package com.github.zimmerlab.gtfcompare.utils;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.model.MappingResult;
import com.github.zimmerlab.gtfcompare.model.TranscriptPair;
import com.github.zimmerlab.gtfcompare.model.comparison.TranscriptComparisonResult;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jgrapht.alg.matching.MaximumWeightBipartiteMatching;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;


public class OverlappingTranscripts {
    private final static Logger logger = LogManager.getLogger(OverlappingTranscripts.class);

    private enum FeatureType {UTR5, CDS, UTR3, INTRON, START_CODON, STOP_CODON}

    private static final double MIN_OVERLAP_FRACTION = 0.00;
    private static final double ENSEMBLE_ALPHA = 0.5;
    private static final double MIN_ENSEMBLE_SCORE = 0.0;
    private static final double IDENTITY_EDGE_WEIGHT = 2.0;

    public static MappingResult<TranscriptPair, TranscriptFeature> map(GtfFile targetGtfFile, GtfFile queryGtfFile) {

        // Build interval trees
        logger.info("Building Interval Trees");
        Map<String, IntervalTreeMap<List<TranscriptFeature>>> targetTrees = buildIntervalTrees(targetGtfFile);
        Map<String, IntervalTreeMap<List<TranscriptFeature>>> queryTrees = buildIntervalTrees(queryGtfFile);

        // Find overlaps
        logger.info("Finding Overlaps");
        List<TranscriptPair> allPairs = findOverlaps(targetTrees, queryTrees).toList();

        // Cluster loci
        logger.info("Cluster Loci");
        Map<TranscriptFeature, List<TranscriptPair>> clusters = clusterOverlappingPairs(allPairs);

        // Prepare final mapping
        logger.info("Compute Global Best Hits via Bipartite Matching");
        var finalMapping = new ConcurrentHashMap<TranscriptFeature, TranscriptFeature>();

        clusters.values().parallelStream().forEach(locusPairs -> {

            // Collect targets and queries

            Set<TranscriptFeature> targets = locusPairs.stream().map(TranscriptPair::getTargetTranscript).collect(Collectors.toCollection(LinkedHashSet::new));
            Set<TranscriptFeature> queries = locusPairs.stream().map(TranscriptPair::getQueryTranscript).collect(Collectors.toCollection(LinkedHashSet::new));
            // Build unified clique for vector construction
            List<TranscriptFeature> clique = Stream.concat(targets.stream(), queries.stream()).collect(Collectors.toList());
            var vectors = buildModelVectorsForClique(clique);

            // Build bipartite graph
            var graph = new SimpleWeightedGraph<TranscriptFeature, DefaultWeightedEdge>(DefaultWeightedEdge.class);
            targets.forEach(graph::addVertex);
            queries.forEach(graph::addVertex);
            //addIdentityMatchEdges(graph, targets, queries);

            for (TranscriptFeature t : targets) {
                for (TranscriptFeature q : queries) {
                    if (overlapFraction(t, q) < MIN_OVERLAP_FRACTION) continue;
                    double jac = jaccardSimilarity(vectors.get(t), vectors.get(q));
                    double chn = chainSimilarity(t, q);
                    double score = ENSEMBLE_ALPHA * jac + (1 - ENSEMBLE_ALPHA) * chn;
                    if (score < MIN_ENSEMBLE_SCORE) continue;
                    DefaultWeightedEdge edge = graph.addEdge(t, q);
                    if (edge != null) {
                        graph.setEdgeWeight(edge, score);
                    }
                }
            }

            //  maximum-weight bipartite matching
            var mwbm = new MaximumWeightBipartiteMatching<>(graph, targets, queries);
            var matching = mwbm.getMatching().getEdges();

            // integrate matches into final mapping
            var cliqueMapping = new LinkedHashMap<TranscriptFeature, TranscriptFeature>();
            for (DefaultWeightedEdge e : matching) {
                TranscriptFeature src = graph.getEdgeSource(e);
                TranscriptFeature tgt = graph.getEdgeTarget(e);
                cliqueMapping.put(src, tgt);
            }

            var matchedTargets = cliqueMapping.keySet();
            var matchedQueries = new HashSet<>(cliqueMapping.values());

            var unmatchedTargets = targets.stream().filter(t -> !matchedTargets.contains(t)).toList();
            var unmatchedQueries = queries.stream().filter(q -> !matchedQueries.contains(q)).toList();

            if (unmatchedTargets.size() == 1 && unmatchedQueries.size() == 1) {
                TranscriptFeature t0 = unmatchedTargets.get(0);
                TranscriptFeature q0 = unmatchedQueries.get(0);
                cliqueMapping.put(t0, q0);
            }

            logger.debug("Cluster: |T|={} |Q|={} edges={}", targets.size(), queries.size(), graph.edgeSet().size());

            finalMapping.putAll(cliqueMapping);
        });

        // Assemble MappingResult
        var mapping = finalMapping.entrySet().stream().map(e -> new TranscriptPair(e.getKey(), e.getValue(), new TranscriptComparisonResult())).toList();

        var allTargets = getAllTranscripts(targetGtfFile);
        var allQueries = getAllTranscripts(queryGtfFile);
        var mappedTargets = finalMapping.keySet();
        var mappedQueries = new HashSet<>(finalMapping.values());

        var unmappedTargets = allTargets.stream().filter(t -> !mappedTargets.contains(t)).toList();
        var unmappedQueries = allQueries.stream().filter(q -> !mappedQueries.contains(q)).toList();
        logger.info("Finished Mapping");
        return new MappingResult<>(mapping, unmappedTargets, unmappedQueries);
    }

    private static Set<TranscriptFeature> getAllTranscripts(GtfFile gtfFile) {
        return gtfFile.getAllGeneFeatureIds().stream().flatMap(geneId -> gtfFile.getGeneFeature(geneId).getTranscripts().stream()).collect(Collectors.toCollection(LinkedHashSet::new));
    }

    private static void addIdentityMatchEdges(SimpleWeightedGraph<TranscriptFeature, DefaultWeightedEdge> graph, Set<TranscriptFeature> targets, Set<TranscriptFeature> queries) {

        for (TranscriptFeature t : targets) {
            for (TranscriptFeature q : queries) {
                if (t.getTranscriptId().equals(q.getTranscriptId())) {
                    var e = graph.addEdge(t, q);
                    if (e != null) {
                        graph.setEdgeWeight(e, IDENTITY_EDGE_WEIGHT);
                    }
                }
            }
        }
    }


    private static Map<String, IntervalTreeMap<List<TranscriptFeature>>> buildIntervalTrees(GtfFile gtfFile) {
        var trees = new HashMap<String, IntervalTreeMap<List<TranscriptFeature>>>();
        Function<TranscriptFeature, String> keyOf = t -> t.getBaseData().getContig() + "_" + (t.getBaseData().isForwardStrand() ? "p" : "n");
        extractTranscripts(gtfFile, trees, keyOf);
        return trees;
    }

    private static Map<TranscriptFeature, List<TranscriptPair>> clusterOverlappingPairs(List<TranscriptPair> allPairs) {
        var ds = new OverlappingTranscripts.DisjointSet<TranscriptFeature>();
        allPairs.forEach(p -> ds.union(p.getTargetTranscript(), p.getQueryTranscript()));

        return allPairs.stream().collect(Collectors.groupingBy(p -> ds.find(p.getTargetTranscript()), LinkedHashMap::new, Collectors.toList()));
    }

    private static List<Integer> getExonLengthChain(TranscriptFeature tf) {
        var exons = tf.getFeatures().stream().filter(seg -> "exon".equals(seg.getBaseData().getType())).sorted(Comparator.comparingInt(s -> s.getBaseData().getStart())).toList();

        var lengths = new ArrayList<Integer>(exons.size());
        for (GtfFeature exon : exons) {
            int len = exon.getBaseData().getEnd() - exon.getBaseData().getStart() + 1;
            lengths.add(len);
        }
        return lengths;
    }

    // 2) Levenshtein-Edit-Distance zwischen zwei Integer-Listen
    private static int editDistance(List<Integer> a, List<Integer> b) {
        int n = a.size(), m = b.size();
        int[][] dp = new int[n + 1][m + 1];
        for (int i = 0; i <= n; i++) dp[i][0] = i;
        for (int j = 0; j <= m; j++) dp[0][j] = j;

        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                int cost = a.get(i - 1).equals(b.get(j - 1)) ? 0 : 1;
                dp[i][j] = Math.min(Math.min(dp[i - 1][j] + 2,      // deletion
                                dp[i][j - 1] + 2),     // insertion
                        dp[i - 1][j - 1] + cost           // substitution
                );
            }
        }
        return dp[n][m];
    }

    private static double chainSimilarity(TranscriptFeature t, TranscriptFeature q) {
        List<Integer> chainT = getExonLengthChain(t);
        List<Integer> chainQ = getExonLengthChain(q);

        if (chainT.isEmpty() && chainQ.isEmpty()) {
            return 1.0;
        }
        int dist = editDistance(chainT, chainQ);
        int maxLen = Math.max(chainT.size(), chainQ.size());
        return maxLen > 0 ? 1.0 - (double) dist / maxLen : 0.0;
    }

    private static double overlapFraction(TranscriptFeature t, TranscriptFeature q) {
        var tStart = t.getBaseData().getStart();
        var tEnd = t.getBaseData().getEnd();
        var qStart = q.getBaseData().getStart();
        var qEnd = q.getBaseData().getEnd();

        var overlapStart = Math.max(tStart, qStart);
        var overlapEnd = Math.min(tEnd, qEnd);
        if (overlapEnd < overlapStart) {
            return 0.0;
        }
        var overlapLen = overlapEnd - overlapStart + 1;
        var tLen = tEnd - tStart + 1;
        var qLen = qEnd - qStart + 1;
        var maxLen = Math.max(tLen, qLen);
        return (double) overlapLen / maxLen;
    }

    // jac 15920
    // edit dist 16499
    // both 16515
    private static Map<TranscriptFeature, TranscriptFeature> mapReciprocalBestHits(Set<TranscriptFeature> targets, Set<TranscriptFeature> queries, Map<TranscriptFeature, EnumMap<OverlappingTranscripts.FeatureType, BitSet>> vectors) {
        var bestQPerT = new HashMap<TranscriptFeature, TranscriptFeature>();
        var bestScoreQ = new HashMap<TranscriptFeature, Double>();

        for (TranscriptFeature t : targets) {
            bestScoreQ.put(t, -1.0);
            for (TranscriptFeature q : queries) {
                if (overlapFraction(t, q) < MIN_OVERLAP_FRACTION) continue;

                // compute both metrics
                double jac = jaccardSimilarity(vectors.get(t), vectors.get(q));
                double chn = chainSimilarity(t, q);

                // ensemble
                double score = ENSEMBLE_ALPHA * jac + (1 - ENSEMBLE_ALPHA) * chn;
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
            for (TranscriptFeature t : targets) {
                if (overlapFraction(t, q) < MIN_OVERLAP_FRACTION) continue;
                if (q.getTranscriptId().equals("ENST00000850868") && q.getTranscriptId().equals("ENST00000458258")) {
                    var a = 2;
                }

                double jac = jaccardSimilarity(vectors.get(t), vectors.get(q));
                double chn = chainSimilarity(t, q);
                double score = ENSEMBLE_ALPHA * jac + (1 - ENSEMBLE_ALPHA) * chn;

                if (score > bestScoreT.get(q)) {
                    bestScoreT.put(q, score);
                    bestTPerQ.put(q, t);
                }
            }
        }

        // 3) reziproke Hits
        var mapping = new LinkedHashMap<TranscriptFeature, TranscriptFeature>();
        for (var e : bestQPerT.entrySet()) {
            TranscriptFeature t = e.getKey();
            TranscriptFeature q = e.getValue();
            if (t.getTranscriptId().equals("ENST00000850868") || q.getTranscriptId().equals("ENST00000850868")) {
                var a = 2;
            }
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
            // intersection
            BitSet ai = (BitSet) a.get(ft).clone();
            ai.and(b.get(ft));
            inter += ai.cardinality();

            // union
            BitSet au = (BitSet) a.get(ft).clone();
            au.or(b.get(ft));
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

    private record FeatureInterval(int start, int end) {
    }

    private static Map<TranscriptFeature, EnumMap<FeatureType, BitSet>> buildModelVectorsForClique(List<TranscriptFeature> clique) {
        // 1. Collect all boundaries (start and end+1) of all segments in the clique
        var coords = new TreeSet<Integer>();
        for (TranscriptFeature tf : clique) {
            for (GtfFeature seg : tf.getFeatures()) {
                coords.add(seg.getBaseData().getStart());
                coords.add(seg.getBaseData().getEnd() + 1);
            }
        }

        // 2. Build sorted list of boundaries and derive disjoint intervals
        var sorted = new ArrayList<>(coords);
        Collections.sort(sorted);
        List<FeatureInterval> subIntervals = new ArrayList<>();
        for (int i = 0; i < sorted.size() - 1; i++) {
            int s = sorted.get(i);
            int e = sorted.get(i + 1) - 1;
            subIntervals.add(new FeatureInterval(s, e));
        }

        // 3. Prepare result map
        Map<TranscriptFeature, EnumMap<FeatureType, BitSet>> vectors = new LinkedHashMap<>();

        // 4. For each transcript, assign each subinterval to its covering segment
        for (TranscriptFeature tf : clique) {
            // initialize an empty BitSet for each FeatureType
            EnumMap<FeatureType, BitSet> map = new EnumMap<>(FeatureType.class);
            for (FeatureType ft : FeatureType.values()) {
                map.put(ft, new BitSet(subIntervals.size()));
            }

            // for each subinterval, determine which segment covers it (if any)
            for (int idx = 0; idx < subIntervals.size(); idx++) {
                FeatureInterval iv = subIntervals.get(idx);
                // find the first segment in tf that spans this entire interval
                for (GtfFeature seg : tf.getFeatures()) {
                    int segS = seg.getBaseData().getStart();
                    int segE = seg.getBaseData().getEnd();
                    if (segS <= iv.start() && segE >= iv.end()) {
                        FeatureType ft = getFeatureType(seg);
                        if (ft != null) {
                            map.get(ft).set(idx);
                        }
                        break; // only one segment can cover
                    }
                }
            }

            vectors.put(tf, map);
        }

        return vectors;
    }

    private static FeatureType getFeatureType(GtfFeature seg) {
        String type = seg.getBaseData().getType();

        if (Constants.UTR5.equals(type)) {
            return FeatureType.UTR5;
        } else if (Constants.CDS.equals(type)) {
            return FeatureType.CDS;
        } else if (Constants.UTR3.equals(type)) {
            return FeatureType.UTR3;
        } else if (Constants.INTRON.equals(type)) {
            return FeatureType.INTRON;
        } else if (Constants.START_CODON.equals(type)) {
            return FeatureType.START_CODON;
        } else if (Constants.STOP_CODON.equals(type)) {
            return FeatureType.STOP_CODON;
        } else {
            return null;
        }
    }
}
