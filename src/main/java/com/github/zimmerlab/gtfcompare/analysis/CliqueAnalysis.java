package com.github.zimmerlab.gtfcompare.analysis;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.model.DisjointSet;
import com.github.zimmerlab.gtfcompare.model.TranscriptPair;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class CliqueAnalysis {
    private final static Logger logger = LogManager.getLogger(CliqueAnalysis.class);

    public void analyze(GtfFile targetGtfFile, GtfFile queryGtfFile, boolean useStrandInKey) {
        setTranscriptPositions(targetGtfFile);
        setTranscriptPositions(queryGtfFile);

        // Build interval trees
        logger.info("Building Interval Trees");
        Map<String, IntervalTreeMap<List<TranscriptFeature>>> targetTrees = buildIntervalTrees(targetGtfFile, useStrandInKey);
        Map<String, IntervalTreeMap<List<TranscriptFeature>>> queryTrees = buildIntervalTrees(queryGtfFile, useStrandInKey);

        // Find overlaps
        logger.info("Finding Overlaps");
        var allPairs = findOverlaps(targetTrees, queryTrees).toList();

        // Cluster loci
        logger.info("Cluster Loci");
        var clusters = clusterOverlappingPairs(allPairs);

    }

    private static Set<TranscriptFeature> getAllTranscripts(GtfFile gtfFile) {
        return gtfFile.getAllGeneFeatureIds().stream().flatMap(geneId -> gtfFile.getGeneFeature(geneId).getTranscripts().stream()).collect(Collectors.toCollection(LinkedHashSet::new));
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

    private static void setTranscriptPositions(GtfFile gtfFile) {
        var transcripts = getAllTranscripts(gtfFile);
        for(var transcript : transcripts){
            var transcriptBaseData = transcript.getBaseData();

            var start = transcriptBaseData.getStart();

            if(start == -1){
                start = transcript.getFeatures().stream().mapToInt(target -> target.getBaseData().getStart())
                        .min()
                        .orElse(Integer.MAX_VALUE);

                transcriptBaseData.setStart(start);
            }

            var stop =  transcriptBaseData.getEnd();
            if(stop == -1){
                stop = transcript.getFeatures().stream().mapToInt(target -> target.getBaseData().getEnd())
                        .max()
                        .orElse(Integer.MIN_VALUE);

                transcriptBaseData.setEnd(stop);
            }
        }
    }


    private static Map<String, IntervalTreeMap<List<TranscriptFeature>>> buildIntervalTrees(GtfFile gtfFile, boolean useStrandInKey) {
        var trees = new HashMap<String, IntervalTreeMap<List<TranscriptFeature>>>();
        Function<TranscriptFeature, String> keyOf;
        if(useStrandInKey) {
            keyOf = t -> t.getBaseData().getContig() + "_" + (t.getBaseData().isForwardStrand() ? "p" : "n");

        } else {
            keyOf = t -> t.getBaseData().getContig();
        }
        extractTranscripts(gtfFile, trees, keyOf);
        return trees;
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

    private static List<List<TranscriptPair>> clusterOverlappingPairs(List<TranscriptPair> allPairs) {
        var ds = new DisjointSet<TranscriptFeature>();
        allPairs.forEach(p -> ds.union(p.getTargetTranscript(), p.getQueryTranscript()));

        var tmp = new HashMap<TranscriptFeature, List<TranscriptPair>>();
        for (TranscriptPair p : allPairs) {
            var representative = ds.find(p.getTargetTranscript());
            tmp.computeIfAbsent(representative, k -> new ArrayList<>()).add(p);
        }
        return new ArrayList<>(tmp.values());

    }
}
