package com.github.zimmerlab.gtfcompare.analysis.clique;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.GeneFeature;
import com.github.zimmerlab.gtfcompare.model.DisjointSet;
import com.github.zimmerlab.gtfcompare.model.GenePair;
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

public class CliqueAnalysis {
    private final static Logger logger = LogManager.getLogger(CliqueAnalysis.class);

    private int genesWithExactMatch;
    private int overallTargetGeneCount;
    private int overallQueryGeneCount;
    private int overallClusters;
    private int unmatchedTargetGenes, unmatchedQueryGenes;
    private int targetGenesInClustersWithoutAnyExact, queryGenesInClustersWithoutAnyExact;
    private final Set<String> exactMatchedGeneIds = new HashSet<>();
    private final Set<String> genesSeenInClustersTarget = new HashSet<>();
    private final Set<String> genesSeenInClustersQuery = new HashSet<>();


    public AnalysisResult analyze(GtfFile targetGtf, GtfFile queryGtf, boolean useStrandInKey) throws IOException {
        var targetTrees = buildIntervalTrees(targetGtf, useStrandInKey);
        var pairs = findOverlaps(targetTrees, queryGtf, useStrandInKey).toList();
        var clusters = clusterOverlappingPairs(pairs);

        computeStatsFromClusters(clusters, targetGtf, queryGtf);

        printUnmatchedGenes(targetGtf, queryGtf);

        return new AnalysisResult(
                targetGtf.getParsedContig(),
                overallTargetGeneCount,
                overallQueryGeneCount,
                overallClusters,
                genesWithExactMatch,
                unmatchedTargetGenes,
                unmatchedQueryGenes,
                targetGenesInClustersWithoutAnyExact,
                queryGenesInClustersWithoutAnyExact
        );
    }

    private void computeStatsFromClusters(List<List<GenePair>> geneClusters, GtfFile targetGtfFile, GtfFile queryGtfFile) {
        overallTargetGeneCount = targetGtfFile.getAllGeneFeatureIds().size();
        overallQueryGeneCount = queryGtfFile.getAllGeneFeatureIds().size();
        overallClusters = geneClusters.size();

        for (var cluster : geneClusters) {
            var targetIds = cluster.stream().map(p -> p.getTargetGene().getGeneId()).collect(Collectors.toCollection(LinkedHashSet::new));
            var queryIds = cluster.stream().map(p -> p.getQueryGene().getGeneId()).collect(Collectors.toCollection(LinkedHashSet::new));

            genesSeenInClustersTarget.addAll(targetIds);
            genesSeenInClustersQuery.addAll(queryIds);

            var exactHere = new HashSet<>(targetIds);
            exactHere.retainAll(queryIds);
            exactMatchedGeneIds.addAll(exactHere);
        }

        unmatchedTargetGenes = overallTargetGeneCount - genesSeenInClustersTarget.size();
        unmatchedQueryGenes = overallQueryGeneCount - genesSeenInClustersQuery.size();

        targetGenesInClustersWithoutAnyExact = (int) genesSeenInClustersTarget.stream().filter(g -> !exactMatchedGeneIds.contains(g)).count();

        queryGenesInClustersWithoutAnyExact = (int) genesSeenInClustersQuery.stream().filter(g -> !exactMatchedGeneIds.contains(g)).count();

        genesWithExactMatch = exactMatchedGeneIds.size();
    }

    private void printUnmatchedGenes(GtfFile targetGtfFile, GtfFile queryGtfFile) throws IOException {
        StringBuilder unmatchedGenes = new StringBuilder();
        var currentContig =targetGtfFile.getParsedContig();
        var allTarget = new HashSet<>(targetGtfFile.getAllGeneFeatureIds());
        var allQuery  = new HashSet<>(queryGtfFile.getAllGeneFeatureIds());

        var unmatchedTarget = new HashSet<>(allTarget);
        unmatchedTarget.removeAll(genesSeenInClustersTarget);
        for (var id : unmatchedTarget) {
            unmatchedGenes.append(String.format("%s\tTARGET\tNO_CLUSTER\t%s%n", currentContig, id));
        }

        var inClusterNoExactTarget = new HashSet<>(genesSeenInClustersTarget);
        inClusterNoExactTarget.removeAll(exactMatchedGeneIds);
        for (var id : inClusterNoExactTarget) {
            unmatchedGenes.append(String.format("%s\tTARGET\tIN_CLUSTER_NO_EXACT\t%s%n", currentContig,id));
        }

        var unmatchedQuery = new HashSet<>(allQuery);
        unmatchedQuery.removeAll(genesSeenInClustersQuery);
        for (var id : unmatchedQuery) {
            unmatchedGenes.append(String.format("%s\tQUERY\tNO_CLUSTER\t%s%n", currentContig, id));
        }

        var inClusterNoExactQuery = new HashSet<>(genesSeenInClustersQuery);
        inClusterNoExactQuery.removeAll(exactMatchedGeneIds);
        for (var id : inClusterNoExactQuery) {
            unmatchedGenes.append(String.format("%s\tQUERY\tIN_CLUSTER_NO_EXACT\t%s%n", currentContig,id));
        }

        Files.writeString(
                Paths.get("output/cliqueAnalysis/unmatched.tsv"),
                unmatchedGenes.toString(),
                StandardCharsets.UTF_8,
                StandardOpenOption.CREATE,
                StandardOpenOption.APPEND
        );    }

    private static Set<GeneFeature> getAllGenes(GtfFile gtfFile) {
        return gtfFile.getAllGeneFeatureIds().stream().map(gtfFile::getGeneFeature).collect(Collectors.toSet());
    }

    private static Stream<GenePair> findOverlaps(Map<String, IntervalTreeMap<List<GeneFeature>>> targetTrees, GtfFile queryGtfFile, boolean useStrandInKey) {
        Function<GeneFeature, String> keyOf = useStrandInKey ? (t -> t.getBaseData().getContig() + "_" + (t.getBaseData().isForwardStrand() ? "p" : "n")) : (t -> t.getBaseData().getContig());

        return getAllGenes(queryGtfFile).stream().flatMap(q -> {
            String key = keyOf.apply(q);
            var tTree = targetTrees.get(key);
            if (tTree == null) return Stream.empty();

            var bd = q.getBaseData();
            var qIv = new Interval(bd.getContig(), bd.getStart(), bd.getEnd());

            return tTree.getOverlapping(qIv).stream().flatMap(List::stream).map(t -> new GenePair(t, q));
        });
    }

    private static Map<String, IntervalTreeMap<List<GeneFeature>>> buildIntervalTrees(GtfFile gtfFile, boolean useStrandInKey) {
        var trees = new HashMap<String, IntervalTreeMap<List<GeneFeature>>>();
        Function<GeneFeature, String> keyOf;
        if (useStrandInKey) {
            keyOf = t -> t.getBaseData().getContig() + "_" + (t.getBaseData().isForwardStrand() ? "p" : "n");

        } else {
            keyOf = t -> t.getBaseData().getContig();
        }
        extractTranscripts(gtfFile, trees, keyOf);
        return trees;
    }

    private static void extractTranscripts(GtfFile gtfFile, Map<String, IntervalTreeMap<List<GeneFeature>>> trees, Function<GeneFeature, String> keyOf) {
        for (var geneId : gtfFile.getAllGeneFeatureIds()) {
            var gene = gtfFile.getGeneFeature(geneId);
            String key = keyOf.apply(gene);
            var bd = gene.getBaseData();
            var iv = new Interval(bd.getContig(), bd.getStart(), bd.getEnd());

            trees.computeIfAbsent(key, k -> new IntervalTreeMap<>()).compute(iv, (interval, list) -> {
                if (list == null) list = new ArrayList<>();
                list.add(gene);
                return list;
            });

        }
    }

    private static List<List<GenePair>> clusterOverlappingPairs(List<GenePair> allPairs) {
        var ds = new DisjointSet<GeneFeature>();
        allPairs.forEach(p -> ds.union(p.getTargetGene(), p.getQueryGene()));

        var tmp = new HashMap<GeneFeature, List<GenePair>>();
        for (var p : allPairs) {
            var representative = ds.find(p.getTargetGene());
            tmp.computeIfAbsent(representative, k -> new ArrayList<>()).add(p);
        }
        return new ArrayList<>(tmp.values());
    }
}
