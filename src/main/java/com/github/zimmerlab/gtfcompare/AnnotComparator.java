package com.github.zimmerlab.gtfcompare;


import com.github.kleinsamuel.gtfutils.GtfConfig;
import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.GeneFeature;
import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.compare.ComparisonConfig;
import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.ComparisonFeature;
import com.github.zimmerlab.gtfcompare.model.FeaturePair;
import com.github.zimmerlab.gtfcompare.model.GenePair;
import com.github.zimmerlab.gtfcompare.model.TranscriptPair;
import com.github.zimmerlab.gtfcompare.model.comparison.ComparisonResult;
import com.github.zimmerlab.gtfcompare.model.comparison.FeatureComparisonResult;
import com.github.zimmerlab.gtfcompare.model.comparison.RegionComparisonResult;
import com.github.zimmerlab.gtfcompare.model.comparison.TranscriptComparisonResult;
import com.github.zimmerlab.gtfcompare.utils.Constants;
import com.github.zimmerlab.gtfcompare.utils.GenomeSequenceExtractor;
import org.springframework.util.StopWatch;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class AnnotComparator {
    private static final List<StopWatch> stopWatches = Constants.STOP_WATCHES;
    private static final ServiceLoader<ComparisonFeature> loader = ServiceLoader.load(ComparisonFeature.class);
    private final GtfFile targetGtf;
    private final GtfFile queryGtf;
    private final GenomeSequenceExtractor targetSequenceExtractor;
    private final GenomeSequenceExtractor querySequenceExtractor;
    private final ComparisonConfig config;
    private final String outputPath;
    private final List<ComparisonResult> comparisonResults = new ArrayList<>();

    public AnnotComparator(GtfFile targetGtf, GtfFile queryGtf, GenomeSequenceExtractor targetSequenceExtractor, GenomeSequenceExtractor querySequenceExtractor, ComparisonConfig config, String outputPath) {
        this.targetGtf = targetGtf;
        this.queryGtf = queryGtf;
        this.targetSequenceExtractor = targetSequenceExtractor;
        this.querySequenceExtractor = querySequenceExtractor;
        this.config = config;
        this.outputPath = outputPath;
    }

    public void compare() {
        var genePairs = getGenePairs();

        for (var pair : genePairs) {
            var result = new ComparisonResult();
            result.setTargetGeneId(pair.getTargetGene() != null ? pair.getTargetGene().getGeneId() : "");
            result.setQueryGeneId(pair.getQueryGene() != null ? pair.getQueryGene().getGeneId() : "");
            compareGene(pair.getTargetGene(), pair.getQueryGene(), result);
            comparisonResults.add(result);
        }

        writeComparisonResult(comparisonResults);
    }

    private void writeComparisonResult(List<ComparisonResult> comparisonResults) {
        try (var writer = new BufferedWriter(new FileWriter((outputPath)))) {
            writer.write("targetGeneId\tqueryGeneId\tcategory\tdifference\n");
            for (var comparisonResult : comparisonResults) {

                var targetGeneId = comparisonResult.getTargetGeneId();
                var queryGeneId = comparisonResult.getQueryGeneId();

                var geneComparison = comparisonResult.getGeneComparison();

                if (geneComparison.isMissingInQueryFile()) {
                    writer.write(targetGeneId + "\t" + queryGeneId + "\tgene\tmissingInQueryFile\n");
                    continue;
                }



                if (comparisonResult.areSameGene()) {
                    writer.write(targetGeneId + "\t" + queryGeneId + "\n");
                } else {

                    if (geneComparison.isStartDifferent()) {
                        writer.write(targetGeneId + "\t" + queryGeneId + "\tgene\tstart\n");
                    }
                    if (geneComparison.isStopDifferent()) {
                        writer.write(targetGeneId + "\t" + queryGeneId + "\tgene\tstop\n");
                    }
                    if (geneComparison.isStrandDifferent()) {
                        writer.write(targetGeneId + "\t" + queryGeneId + "\tgene\tstrand\n");
                    }
                    if (!geneComparison.getSequenceComparison().isSameSequence()) {
                        writer.write(targetGeneId + "\t" + queryGeneId+ "\tgene\tseq\n");
                    }
                    if (geneComparison.isDifferentLength()) {
                        writer.write(targetGeneId + "\t" + queryGeneId+ "\tgene\tlength\n");
                    }

                    var transcriptComparisonResult = comparisonResult.getTranscriptComparisons();

                    for (var transcriptComparison : transcriptComparisonResult) {
                        if (transcriptComparison.isStartDifferent()) {
                            writer.write(targetGeneId + "\t" + queryGeneId+ "\ttranscript\tstart_" + transcriptComparison.getQueryTranscriptId() + "\n");
                        }
                        if (transcriptComparison.isStopDifferent()) {
                            writer.write(targetGeneId + "\t" + queryGeneId+ "\ttranscript\tstop_" + transcriptComparison.getQueryTranscriptId() + "\n");
                        }
                        if (transcriptComparison.isSequenceDifferent()) {
                            writer.write(targetGeneId + "\t" + queryGeneId+ "\ttranscript\tseq_" + transcriptComparison.getQueryTranscriptId() + "\n");
                        }
                        if (transcriptComparison.isTranscriptMissingInTargetGene()) {
                            writer.write(targetGeneId + "\t" + queryGeneId+ "\ttranscript\tmissingInFile1_" + transcriptComparison.getQueryTranscriptId() + "\n");
                        }
                        if (transcriptComparison.isTranscriptMissingInQueryGene()) {
                            writer.write(targetGeneId + "\t" + queryGeneId+ "\ttranscript\tmissingInFile2_" + transcriptComparison.getQueryTranscriptId() + "\n");
                        }
                        if (transcriptComparison.isLengthDifferent()) {
                            writer.write(targetGeneId + "\t" + queryGeneId+ "\ttranscript\tlength_" + transcriptComparison.getQueryTranscriptId() + "\n");

                        }

                        for (FeatureComparisonResult featureComparison : transcriptComparison.getFeatureComparisons()) {
                            if (featureComparison.isMissingInTargetTranscript()) {
                                writer.write(targetGeneId + "\t" + queryGeneId+ "\tfeature\tmissingInTranscript1_" + transcriptComparison.getQueryTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                            }
                            if (featureComparison.isMissingInQueryTranscript()) {
                                writer.write(targetGeneId + "\t" + queryGeneId+ "\tfeature\tmissingInTranscript2_" + transcriptComparison.getQueryTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                            }

                            for (RegionComparisonResult regionComparison : featureComparison.getRegionComparisons()) {
                                if (regionComparison.isLengthDifferenceFound()) {
                                    writer.write(targetGeneId + "\t" + queryGeneId+ "\tfeature\tlength_" + transcriptComparison.getQueryTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                }
                                if (regionComparison.isPositionDifferenceFound()) {
                                    writer.write(targetGeneId + "\t" + queryGeneId+ "\tfeature\tposition_" + transcriptComparison.getQueryTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                }
                                if (regionComparison.isSequenceDifferenceFound()) {
                                    writer.write(targetGeneId + "\t" + queryGeneId+ "\tfeature\tseq_" + transcriptComparison.getQueryTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                }
                                if (regionComparison.isMissingInTargetFile()) {
                                    writer.write(targetGeneId + "\t" + queryGeneId+ "\tfeature\tmissingFeatureEntryFile1_" + transcriptComparison.getQueryTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                }
                                if (regionComparison.isMissingInQueryFile()) {
                                    writer.write(targetGeneId + "\t" + queryGeneId+ "\tfeature\tmissingFeatureEntryFile2_" + transcriptComparison.getQueryTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                }
                            }
                        }
                    }
                }
            }

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private List<GenePair> getGenePairs() {
        var targetGeneMap = new HashMap<String, GeneFeature>();
        var queryGeneMap = new HashMap<String, GeneFeature>();

        for (var geneId : targetGtf.getAllGeneFeatureIds()) {
            targetGeneMap.put(geneId, targetGtf.getGeneFeature(geneId));
        }

        for (var geneId : queryGtf.getAllGeneFeatureIds()) {
            queryGeneMap.put(geneId, queryGtf.getGeneFeature(geneId));
        }

        var allGeneIds = new HashSet<String>();
        allGeneIds.addAll(targetGeneMap.keySet());
        allGeneIds.addAll(queryGeneMap.keySet());

        var genePairs = new ArrayList<GenePair>();
        for (var geneId : allGeneIds) {
            GeneFeature t1 = targetGeneMap.get(geneId);
            GeneFeature t2 = queryGeneMap.get(geneId);
            genePairs.add(new GenePair(t1, t2));
        }
        return genePairs;
    }

    private void compareGene(GeneFeature targetGene, GeneFeature queryGene, ComparisonResult result) {

        if (targetGene == null) {
            result.setAreSameGene(false);
            result.getGeneComparison().setMissingInTargetFile(true);
            return;
        }
        if (queryGene == null) {
            result.setAreSameGene(false);
            result.getGeneComparison().setMissingInQueryFile(true);
            return;
        }

        var transcriptPairs = getTranscriptPairs(targetGene, queryGene, result);
        for (var tp : transcriptPairs) {
            compareTranscript(tp, result);
        }
    }

    private void compareTranscript(TranscriptPair tp, ComparisonResult geneResult) {

        var transcriptComparisonResult = tp.getTranscriptComparisonResult();
        var targetTranscript = tp.getTargetTranscript();
        var queryTranscript = tp.getQueryTranscript();

        if (targetTranscript == null || queryTranscript == null) {
            handleMissingTranscript(targetTranscript, queryTranscript, transcriptComparisonResult, geneResult);
            return;
        }

        transcriptComparisonResult.setTargetTranscriptId(targetTranscript.getTranscriptId());
        transcriptComparisonResult.setQueryTranscriptId(queryTranscript.getTranscriptId());

        var targetMap = mapFeaturesByType(targetTranscript);
        var queryMap = mapFeaturesByType(queryTranscript);

        compareFeatures(targetMap, queryMap, transcriptComparisonResult, geneResult);

    }

    private void handleMissingTranscript(TranscriptFeature a, TranscriptFeature b, TranscriptComparisonResult txResult, ComparisonResult geneResult) {

        geneResult.setAreSameGene(false);
        txResult.setAreSameTranscript(false);

        if (a == null) {
            txResult.setTranscriptMissingInTargetGene(true);
            txResult.setQueryTranscriptId(b.getTranscriptId());
        } else {
            txResult.setTranscriptMissingInQueryGene(true);
            txResult.setTargetTranscriptId(a.getTranscriptId());
        }
    }

    private Map<String, List<GtfFeature>> mapFeaturesByType(TranscriptFeature transcriptFeature) {
        var map = new HashMap<String, List<GtfFeature>>();
        for (var f : transcriptFeature.getFeatures()) {
            var type = GtfConfig.getDefault(f.getBaseData().getType());
            map.computeIfAbsent(type, k -> new ArrayList<>()).add(f);
        }
        return map;
    }

    private void compareFeatures(Map<String, List<GtfFeature>> targetMap, Map<String, List<GtfFeature>> queryMap, TranscriptComparisonResult transcriptComparisonResult, ComparisonResult geneResult) {

        for (var entry : targetMap.entrySet()) {
            var featureType = entry.getKey();
            List<GtfFeature> targets = entry.getValue();
            List<GtfFeature> queries = queryMap.getOrDefault(featureType, List.of());

            var featureComparisonResult = new FeatureComparisonResult();
            featureComparisonResult.setFeatureType(featureType);
            transcriptComparisonResult.addFeatureComparison(featureComparisonResult);

            if (queries.isEmpty()) {
                featureComparisonResult.setAreSameFeatures(false);
                featureComparisonResult.setMissingInQueryTranscript(true);
                transcriptComparisonResult.setAreSameTranscript(false);
                geneResult.setAreSameGene(false);
                continue;
            }

            List<FeaturePair> pairs = pairByPosition(targets, queries);
            for (var pair : pairs) {
                compareFeaturePair(pair, featureComparisonResult, transcriptComparisonResult, geneResult);
            }
        }

        for (var featureName : queryMap.keySet()) {
            if (!targetMap.containsKey(featureName)) {
                var featRes = new FeatureComparisonResult();
                featRes.setFeatureType(featureName);
                featRes.setAreSameFeatures(false);
                transcriptComparisonResult.addFeatureComparison(featRes);
                featRes.setMissingInTargetTranscript(true);
                transcriptComparisonResult.setAreSameTranscript(false);
                geneResult.setAreSameGene(false);
            }
        }
    }

    private void compareFeaturePair(FeaturePair pair, FeatureComparisonResult featureComparisonResult, TranscriptComparisonResult transcriptComparisonResult, ComparisonResult geneResult) {

        var targetFeature = pair.getTarget();
        var queryFeature = pair.getQuery();

        if (targetFeature == null || queryFeature == null) {
            markMissingFeature(targetFeature, queryFeature, featureComparisonResult, transcriptComparisonResult, geneResult);
            return;
        }

        var targetBaseData = targetFeature.getBaseData();
        var queryBaseData = queryFeature.getBaseData();

        var regionComparisonResult = new RegionComparisonResult(targetBaseData.getStart(), targetBaseData.getEnd(), queryBaseData.getStart(), queryBaseData.getEnd());
        featureComparisonResult.addRegionComparison(regionComparisonResult);
        var ctx = new ComparisonContext(targetFeature, queryFeature, config, targetSequenceExtractor, querySequenceExtractor);
        for (var comp : loader) {
            if (!config.isEnabled(comp.getName()))
                continue;

            var changed = comp.compare(ctx);

            if (changed) {
                addToRegionComparison(comp.getName(), regionComparisonResult, featureComparisonResult, transcriptComparisonResult, geneResult);
            }

            System.out.printf(
                    "Feature %s: %s changed: %b%n",
                    targetFeature.getBaseData().getType(), comp.getName(), changed
            );
        }
    }

    private void addToRegionComparison(String name, RegionComparisonResult regionComparisonResult, FeatureComparisonResult featureComparisonResult ,TranscriptComparisonResult transcriptComparisonResult, ComparisonResult geneResult) {
        switch (name) {
            case Constants.LENGTH_COMPARATOR_NAME:
                transcriptComparisonResult.setAreSameTranscript(false);
                geneResult.setAreSameGene(false);
                featureComparisonResult.setAreSameFeatures(false);
                regionComparisonResult.setLengthDifferenceFound(true);
                break;
            case Constants.SEQUENCE_COMPARATOR_NAME:
                transcriptComparisonResult.setAreSameTranscript(false);
                geneResult.setAreSameGene(false);
                featureComparisonResult.setAreSameFeatures(false);
                regionComparisonResult.setSequenceDifferenceFound(true);
                // TODO Add Sequence Difference
                break;
            case Constants.START_COMPARATOR_NAME:
                transcriptComparisonResult.setAreSameTranscript(false);
                geneResult.setAreSameGene(false);
                featureComparisonResult.setAreSameFeatures(false);
                regionComparisonResult.setStartDifferent(true);
                break;
            case Constants.STOP_COMPARATOR_NAME:
                transcriptComparisonResult.setAreSameTranscript(false);
                geneResult.setAreSameGene(false);
                featureComparisonResult.setAreSameFeatures(false);
                regionComparisonResult.setEndDifferent(true);
                break;
            default:
                // TODO add logging
                break;
        }
    }

    private void markMissingFeature(
            GtfFeature targetFeature, GtfFeature queryFeature,
            FeatureComparisonResult featRes, TranscriptComparisonResult transcriptComparisonResult,
            ComparisonResult geneResult) {

        geneResult.setAreSameGene(false);
        transcriptComparisonResult.setAreSameTranscript(false);
        featRes.setAreSameFeatures(false);
        if (targetFeature == null) {
            var bd = queryFeature.getBaseData();
            var regionComparison = new RegionComparisonResult(-1, -1, bd.getStart(), bd.getEnd());
            regionComparison.setMissingInTargetFile(true);
            featRes.addRegionComparison(regionComparison);
        } else {
            var bd = targetFeature.getBaseData();
            var regionComparison = new RegionComparisonResult(bd.getStart(), bd.getEnd(), -1, -1);
            regionComparison.setMissingInQueryFile(true);
            featRes.addRegionComparison(regionComparison);
        }
    }


    private List<FeaturePair> pairByPosition(
            List<GtfFeature> targets,
            List<GtfFeature> queries) {

        Comparator<GtfFeature> byStartThenEnd = Comparator
                .comparingInt((GtfFeature f) -> f.getBaseData().getStart())
                .thenComparingInt(f -> f.getBaseData().getEnd());

        var a = new ArrayList<>(targets);
        var b = new ArrayList<>(queries);
        a.sort(byStartThenEnd);
        b.sort(byStartThenEnd);

        var n = Math.max(a.size(), b.size());
        var pairs = new ArrayList<FeaturePair>(n);
        for (int i = 0; i < n; i++) {
            GtfFeature ta = i < a.size() ? a.get(i) : null;
            GtfFeature qb = i < b.size() ? b.get(i) : null;
            pairs.add(new FeaturePair(ta, qb));
        }
        return pairs;
    }

    private List<TranscriptPair> getTranscriptPairs(GeneFeature targetGene, GeneFeature queryGene, ComparisonResult comparisonResult) {
        var stopWatch = new StopWatch();
        stopWatches.add(stopWatch);
        stopWatch.start("getTranscriptPairs");
        var targetTranscriptMap = new HashMap<String, TranscriptFeature>();
        var queryTranscriptMap = new HashMap<String, TranscriptFeature>();

        for (TranscriptFeature transcript : targetGene.getTranscripts()) {
            targetTranscriptMap.put(transcript.getTranscriptId(), transcript);
        }

        for (TranscriptFeature transcript : queryGene.getTranscripts()) {
            queryTranscriptMap.put(transcript.getTranscriptId(), transcript);
        }

        var allTranscriptIds = new HashSet<String>();
        allTranscriptIds.addAll(targetTranscriptMap.keySet());
        allTranscriptIds.addAll(queryTranscriptMap.keySet());

        var transcriptPairs = new ArrayList<TranscriptPair>();
        for (var transcriptId : allTranscriptIds) {
            var transcriptComparisonResult = new TranscriptComparisonResult();
            comparisonResult.addTranscriptComparison(transcriptComparisonResult);
            TranscriptFeature t1 = targetTranscriptMap.get(transcriptId);
            TranscriptFeature t2 = queryTranscriptMap.get(transcriptId);
            transcriptPairs.add(new TranscriptPair(t1, t2, transcriptComparisonResult));
        }

        stopWatch.stop();
        return transcriptPairs;
    }

}
