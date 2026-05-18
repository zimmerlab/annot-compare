package com.github.zimmerlab.gtfcompare.analysis.clique;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.model.DisjointSet;
import com.github.zimmerlab.gtfcompare.model.TranscriptPair;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class CliqueAnalysisTranscripts {
    private final static Logger logger = LogManager.getLogger(CliqueAnalysisTranscripts.class);

    private int transcriptsWithExactMatch;
    private int overallTargetTranscriptCount;
    private int overallQueryTranscriptCount;
    private int overallClusters;
    private int unmatchedTargetTranscripts, unmatchedQueryTranscripts;
    private int targetTranscriptsInClustersWithoutAnyExact, queryTranscriptsInClustersWithoutAnyExact;
    private final Set<String> exactMatchedTranscriptsIds = new HashSet<>();
    private final Set<String> transcriptsSeenInClustersTarget = new HashSet<>();
    private final Set<String> transcriptsSeenInClustersQuery = new HashSet<>();

    public AnalysisResult analyze(GtfFile targetGtf, GtfFile queryGtf, boolean useStrandInKey, String outputPath) throws IOException {
        var targetTrees = buildIntervalTrees(targetGtf, useStrandInKey);
        var pairs = findOverlaps(targetTrees, queryGtf, useStrandInKey).toList();
        var clusters = clusterOverlappingPairs(pairs);

        computeStatsFromClusters(clusters, targetGtf, queryGtf);

        printUnmatchedGenes(targetGtf, queryGtf, outputPath);

        return new AnalysisResult(
                targetGtf.getParsedContig(),
                overallTargetTranscriptCount,
                overallQueryTranscriptCount,
                overallClusters,
                transcriptsWithExactMatch,
                unmatchedTargetTranscripts,
                unmatchedQueryTranscripts,
                targetTranscriptsInClustersWithoutAnyExact,
                queryTranscriptsInClustersWithoutAnyExact
        );
    }

    private void computeStatsFromClusters(List<List<TranscriptPair>> geneClusters, GtfFile targetGtfFile, GtfFile queryGtfFile) {
        overallTargetTranscriptCount = targetGtfFile.getAllGeneFeatureIds().size();
        overallQueryTranscriptCount = queryGtfFile.getAllGeneFeatureIds().size();
        overallClusters = geneClusters.size();

        for (var cluster : geneClusters) {
            var targetIds = cluster.stream().map(p -> p.getTarget().getTranscriptId()).collect(Collectors.toCollection(LinkedHashSet::new));
            var queryIds = cluster.stream().map(p -> p.getQuery().getTranscriptId()).collect(Collectors.toCollection(LinkedHashSet::new));

            transcriptsSeenInClustersTarget.addAll(targetIds);
            transcriptsSeenInClustersQuery.addAll(queryIds);

            var exactHere = new HashSet<>(targetIds);
            exactHere.retainAll(queryIds);
            exactMatchedTranscriptsIds.addAll(exactHere);
        }

        unmatchedTargetTranscripts = overallTargetTranscriptCount - transcriptsSeenInClustersTarget.size();
        unmatchedQueryTranscripts = overallQueryTranscriptCount - transcriptsSeenInClustersQuery.size();

        targetTranscriptsInClustersWithoutAnyExact = (int) transcriptsSeenInClustersTarget.stream().filter(g -> !exactMatchedTranscriptsIds.contains(g)).count();

        queryTranscriptsInClustersWithoutAnyExact = (int) transcriptsSeenInClustersQuery.stream().filter(g -> !exactMatchedTranscriptsIds.contains(g)).count();

        transcriptsWithExactMatch = exactMatchedTranscriptsIds.size();
    }

    private void printUnmatchedGenes(GtfFile targetGtfFile, GtfFile queryGtfFile, String outputPath) throws IOException {
        StringBuilder unmatchedGenes = new StringBuilder();
        var currentContig =targetGtfFile.getParsedContig();
        var allTarget = new HashSet<>(targetGtfFile.getAllGeneFeatureIds());
        var allQuery  = new HashSet<>(queryGtfFile.getAllGeneFeatureIds());

        var unmatchedTarget = new HashSet<>(allTarget);
        unmatchedTarget.removeAll(transcriptsSeenInClustersTarget);
        for (var id : unmatchedTarget) {
            unmatchedGenes.append(String.format("%s\tTARGET\tNO_CLUSTER\t%s%n", currentContig, id));
        }

        var inClusterNoExactTarget = new HashSet<>(transcriptsSeenInClustersTarget);
        inClusterNoExactTarget.removeAll(exactMatchedTranscriptsIds);
        for (var id : inClusterNoExactTarget) {
            unmatchedGenes.append(String.format("%s\tTARGET\tIN_CLUSTER_NO_EXACT\t%s%n", currentContig,id));
        }

        var unmatchedQuery = new HashSet<>(allQuery);
        unmatchedQuery.removeAll(transcriptsSeenInClustersQuery);
        for (var id : unmatchedQuery) {
            unmatchedGenes.append(String.format("%s\tQUERY\tNO_CLUSTER\t%s%n", currentContig, id));
        }

        var inClusterNoExactQuery = new HashSet<>(transcriptsSeenInClustersQuery);
        inClusterNoExactQuery.removeAll(exactMatchedTranscriptsIds);
        for (var id : inClusterNoExactQuery) {
            unmatchedGenes.append(String.format("%s\tQUERY\tIN_CLUSTER_NO_EXACT\t%s%n", currentContig,id));
        }

        Files.writeString(
                Paths.get(outputPath),
                unmatchedGenes.toString(),
                StandardCharsets.UTF_8,
                StandardOpenOption.CREATE,
                StandardOpenOption.APPEND
        );    }

    private static Stream<TranscriptPair> findOverlaps(Map<String, IntervalTreeMap<List<TranscriptFeature>>> targetTrees, GtfFile queryGtfFile, boolean useStrandInKey) {
        Function<TranscriptFeature, String> keyOf = useStrandInKey ? (t -> t.getBaseData().getContig() + "_" + (t.getBaseData().isForwardStrand() ? "p" : "n")) : (t -> t.getBaseData().getContig());

        var queryTranscripts = queryGtfFile.getAllGeneFeatureIds().stream().flatMap(g -> queryGtfFile.getGeneFeature(g).getTranscripts().stream());

        return queryTranscripts.flatMap(q -> {
            String key = keyOf.apply(q);
            var tTree = targetTrees.get(key);
            if (tTree == null) return Stream.empty();

            var qBaseData = q.getBaseData();
            var qInterval = new Interval(qBaseData.getContig(), qBaseData.getStart(), qBaseData.getEnd());
            return tTree.getOverlapping(qInterval)
                    .stream()
                    .flatMap(List::stream)
                    .map(t -> new TranscriptPair(t, q));
        });
    }

    private static Map<String, IntervalTreeMap<List<TranscriptFeature>>> buildIntervalTrees(GtfFile gtfFile, boolean useStrandInKey) {
        var trees = new HashMap<String, IntervalTreeMap<List<TranscriptFeature>>>();
        Function<TranscriptFeature, String> keyOf;
        if (useStrandInKey) {
            keyOf = t -> t.getBaseData().getContig() + "_" + (t.getBaseData().isForwardStrand() ? "p" : "n");

        } else {
            keyOf = t -> t.getBaseData().getContig();
        }
        extractTranscripts(gtfFile, trees, keyOf);
        return trees;
    }

    private static void extractTranscripts(GtfFile gtfFile, Map<String, IntervalTreeMap<List<TranscriptFeature>>> trees, Function<TranscriptFeature, String> keyOf) {
        for (var geneId : gtfFile.getAllGeneFeatureIds()) {
            var gene = gtfFile.getGeneFeature(geneId);
            for (var transcript : gene.getTranscripts()) {
                String key = keyOf.apply(transcript);
                var bd = transcript.getBaseData();
                var iv = new Interval(bd.getContig(), bd.getStart(), bd.getEnd());

                trees.computeIfAbsent(key, k -> new IntervalTreeMap<>()).compute(iv, (interval, list) -> {
                    if (list == null) list = new ArrayList<>();
                    list.add(transcript);
                    return list;
                });
            }
        }
    }

    private static List<List<TranscriptPair>> clusterOverlappingPairs(List<TranscriptPair> allPairs) {
        var ds = new DisjointSet<TranscriptFeature>();
        allPairs.forEach(p -> ds.union(p.getTarget(), p.getQuery()));

        var tmp = new HashMap<TranscriptFeature, List<TranscriptPair>>();
        for (var p : allPairs) {
            var representative = ds.find(p.getTarget());
            tmp.computeIfAbsent(representative, k -> new ArrayList<>()).add(p);
        }
        return new ArrayList<>(tmp.values());
    }
}
