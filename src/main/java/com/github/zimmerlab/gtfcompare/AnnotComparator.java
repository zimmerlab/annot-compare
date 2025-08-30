package com.github.zimmerlab.gtfcompare;


import com.github.kleinsamuel.gtfutils.GtfConfig;
import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.GeneFeature;
import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.compare.*;
import com.github.zimmerlab.gtfcompare.mapping.ExonMapping;
import com.github.zimmerlab.gtfcompare.model.FeaturePair;
import com.github.zimmerlab.gtfcompare.model.GenePair;
import com.github.zimmerlab.gtfcompare.model.TranscriptPair;
import com.github.zimmerlab.gtfcompare.model.comparison.*;
import com.github.zimmerlab.gtfcompare.utils.Constants;
import com.github.zimmerlab.gtfcompare.utils.GenomeSequenceExtractor;
import com.github.zimmerlab.gtfcompare.utils.ResultWriter;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.util.StopWatch;

import java.util.*;
import java.util.function.BiPredicate;
import java.util.function.ToIntFunction;

import static com.github.zimmerlab.gtfcompare.mapping.ExonMapping.*;

public class AnnotComparator {
    private static final List<StopWatch> stopWatches = Constants.STOP_WATCHES;
    private static final Logger logger = LogManager.getLogger(AnnotComparator.class);
    private static final ServiceLoader<ComparisonFeature> featureLoader = ServiceLoader.load(ComparisonFeature.class);
    private static final ServiceLoader<CDSComparisonFeature> cdsLoader = ServiceLoader.load(CDSComparisonFeature.class);
    private static final ServiceLoader<TranscriptComparisonFeature> transcriptLoader = ServiceLoader.load(TranscriptComparisonFeature.class);
    private static final ServiceLoader<GeneComparisonFeature> geneLoader = ServiceLoader.load(GeneComparisonFeature.class);
    private final GtfFile targetGtf;
    private final GtfFile queryGtf;
    private final GenomeSequenceExtractor targetSequenceExtractor;
    private final GenomeSequenceExtractor querySequenceExtractor;
    private final ComparisonConfig config;
    private final String outputPath;
    private final List<ComparisonResult> comparisonResults = new ArrayList<>();
    private final List<TranscriptPair> transcriptPairs;

    public AnnotComparator(GtfFile targetGtf, GtfFile queryGtf, GenomeSequenceExtractor targetSequenceExtractor, GenomeSequenceExtractor querySequenceExtractor, ComparisonConfig config, String outputPath, List<TranscriptPair> transcriptPairs) {
        this.targetGtf = targetGtf;
        this.queryGtf = queryGtf;
        this.targetSequenceExtractor = targetSequenceExtractor;
        this.querySequenceExtractor = querySequenceExtractor;
        this.config = config;
        this.outputPath = outputPath;
        this.transcriptPairs = transcriptPairs;
    }

    public void compare() {
        //var genePairs = getGenePairs();
        //logger.info("number of gene pairs with new mapping: {}", genePairs.size());

        var overallTranscriptPairs = transcriptPairs.size();

        int idx = 0;
        var everyBiotypeAllowed = config.getAllowedGeneBiotypes().isEmpty();
        for (var pair : transcriptPairs) {
            var isTranscriptBiotypeAllowed = isBiotypeAllowed(pair);
            idx++;
            if (!everyBiotypeAllowed && !isTranscriptBiotypeAllowed.isAllowed) continue;

            var geneBiotypeResponse = getGeneBiotype(pair);
            var result = new ComparisonResult();
            result.getGeneComparison().setQueryBiotype(geneBiotypeResponse.queryBiotype);
            result.getGeneComparison().setTargetBiotype(geneBiotypeResponse.targetBiotype);
            var transcriptComparisonResult = pair.getTranscriptComparisonResult();
            transcriptComparisonResult.setQueryBiotype(isTranscriptBiotypeAllowed.queryBiotype);
            transcriptComparisonResult.setTargetBiotype(isTranscriptBiotypeAllowed.targetBiotype);
            result.addTranscriptComparison(transcriptComparisonResult);
            var queryGeneId = pair.getQueryTranscript().getBaseData().getAttributes("gene_id");
            var targetGeneId = pair.getTargetTranscript().getBaseData().getAttributes("gene_id");
            result.setQueryGeneId(queryGeneId != null ? queryGeneId.getFirst() : "");
            result.setTargetGeneId(targetGeneId != null ? targetGeneId.getFirst() : "");
            compareTranscript(pair, result);

            comparisonResults.add(result);

            if (idx % 100 == 0) {
                logger.info(String.format("%d transcript pairs of %d analyzed", idx, overallTranscriptPairs));
            }
        }

        ResultWriter.writeComparisonResult(comparisonResults, outputPath, config);
    }

    private record BiotypeAllowedResponse(boolean isAllowed, String targetBiotype, String queryBiotype) {
    }

    private BiotypeAllowedResponse isBiotypeAllowed(TranscriptPair pair) {
        var allowedBiotypes = config.getAllowedGeneBiotypes();
        var targetBaseData = pair.getTargetTranscript().getBaseData();
        var queryBaseData = pair.getQueryTranscript().getBaseData();

        var hasTargetBiotypeAttribute = targetBaseData.getAttributes("transcript_biotype") != null;
        var hasQueryBiotypeAttribute = queryBaseData.getAttributes("transcript_biotype") != null;

        var targetBiotype = hasTargetBiotypeAttribute ? targetBaseData.getAttributes("transcript_biotype").getFirst() : targetBaseData.getSource();
        var queryBiotype = hasQueryBiotypeAttribute ? queryBaseData.getAttributes("transcript_biotype").getFirst() : targetBaseData.getSource();

        return new BiotypeAllowedResponse(allowedBiotypes.contains(queryBiotype) || allowedBiotypes.contains(targetBiotype), targetBiotype, queryBiotype);
    }

    private record GetGeneBiotypeResponse(String targetBiotype, String queryBiotype) {
    }

    private GetGeneBiotypeResponse getGeneBiotype(TranscriptPair pair) {
        var targetBaseData = pair.getTargetTranscript().getBaseData();
        var queryBaseData = pair.getQueryTranscript().getBaseData();

        var hasTargetBiotypeAttribute = targetBaseData.getAttributes("gene_biotype") != null;
        var hasQueryBiotypeAttribute = queryBaseData.getAttributes("gene_biotype") != null;

        var targetBiotype = hasTargetBiotypeAttribute ? targetBaseData.getAttributes("gene_biotype").getFirst() : targetBaseData.getSource();
        var queryBiotype = hasQueryBiotypeAttribute ? queryBaseData.getAttributes("gene_biotype").getFirst() : targetBaseData.getSource();

        return new GetGeneBiotypeResponse(targetBiotype, queryBiotype);
    }

    private List<GenePair> getGenePairsByExactId() {
        var targetGeneMap = new HashMap<String, GeneFeature>();
        var queryGeneMap = new HashMap<String, GeneFeature>();

        var everyBiotypeAllowed = config.getAllowedGeneBiotypes().isEmpty();
        for (var geneId : targetGtf.getAllGeneFeatureIds()) {
            var gene = targetGtf.getGeneFeature(geneId);

            if (everyBiotypeAllowed) {
                targetGeneMap.put(geneId, gene);
                continue;
            }
            var baseData = gene.getBaseData();
            var geneBiotypeEmpty = baseData.getAttributes("gene_biotype").isEmpty();


            String geneBiotype = geneBiotypeEmpty ? baseData.getSource() : baseData.getAttributes("gene_biotype").getFirst();
            if (config.getAllowedGeneBiotypes().contains(geneBiotype)) targetGeneMap.put(geneId, gene);
        }

        for (var geneId : queryGtf.getAllGeneFeatureIds()) {
            var gene = queryGtf.getGeneFeature(geneId);
            if (everyBiotypeAllowed) {
                queryGeneMap.put(geneId, gene);
                continue;
            }
            var baseData = gene.getBaseData();
            var geneBiotypeEmpty = baseData.getAttributes("gene_biotype").isEmpty();

            String geneBiotype = geneBiotypeEmpty ? baseData.getSource() : baseData.getAttributes("gene_biotype").getFirst();
            if (config.getAllowedGeneBiotypes().contains(geneBiotype)) queryGeneMap.put(geneId, gene);
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
        if (targetGene == null || queryGene == null) {
            handleMissingGene(targetGene == null, result);
            return;
        }

        var transcriptPairs = getTranscriptPairs(targetGene, queryGene, result);
        for (var tp : transcriptPairs) {
            compareTranscript(tp, result);
        }

        var targetTranscripts = targetGene.getTranscripts();
        var queryTranscripts = queryGene.getTranscripts();

        var ctx = new ComparisonContext(targetGene, queryGene, "", "", config, targetSequenceExtractor, querySequenceExtractor, targetTranscripts, queryTranscripts);

        for (var comp : geneLoader) {
            if (!config.isEnabled(comp.getName())) continue;

            var changed = comp.compare(ctx);

            if (changed) {
                addToGeneComparison(comp.getName(), result);
            }

            logger.debug("Gene {}: {} changed: {}", targetGene.getBaseData().getType(), comp.getName(), changed);
        }

        var targetGeneBiotypeEmpty = targetGene.getBaseData().getAttributes("gene_biotype").isEmpty();
        var queryGeneBiotypeEmpty = queryGene.getBaseData().getAttributes("gene_biotype").isEmpty();
        var geneComparison = result.getGeneComparison();
        if (!targetGeneBiotypeEmpty) {
            geneComparison.setTargetBiotype(targetGene.getBaseData().getAttributes("gene_biotype").getFirst());
        }

        if (!queryGeneBiotypeEmpty) {
            geneComparison.setQueryBiotype(queryGene.getBaseData().getAttributes("gene_biotype").getFirst());
        }

    }

    private void handleMissingGene(boolean isTargetNull, ComparisonResult result) {
        if (isTargetNull) {
            result.setAreSameGene(false);
            var geneComparison = result.getGeneComparison();
            geneComparison.setAreSameGene(false);
            geneComparison.setMissingInTargetFile(true);
        } else {
            result.setAreSameGene(false);
            var geneComparison = result.getGeneComparison();
            geneComparison.setAreSameGene(false);
            geneComparison.setMissingInQueryFile(true);
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

        var pairedExons = ExonMapping.pairExonsByGapAlignment(targetTranscript, queryTranscript, 2, 200, 200, 0.2);

        var targetMap = mapFeaturesByType(targetTranscript);
        var queryMap = mapFeaturesByType(queryTranscript);

        compareFeatures(pairedExons, targetTranscript, queryTranscript, transcriptComparisonResult, geneResult, 10);

        // todo mby checken ob aktiviert in config?
        var targetFeatures = targetTranscript.getFeatures();
        var queryFeatures = queryTranscript.getFeatures();

        var targetBiotype = transcriptComparisonResult.getTargetBiotype();
        var queryBiotype = transcriptComparisonResult.getQueryBiotype();

        var ctx = new ComparisonContext(targetTranscript, queryTranscript, targetBiotype, queryBiotype, config, targetSequenceExtractor, querySequenceExtractor, targetFeatures, queryFeatures);
        for (var comp : transcriptLoader) {
            if (!config.isEnabled(comp.getName())) continue;

            var changed = comp.compare(ctx);

            if (changed) {
                addToTranscriptComparison(comp.getName(), transcriptComparisonResult, geneResult);
            }

            logger.debug("Transcript {}: {} changed: {}", targetTranscript.getBaseData().getType(), comp.getName(), changed);
        }

        transcriptComparisonResult.setQueryStart(ctx.getQueryTranscriptFeaturesMin());
        transcriptComparisonResult.setQueryStop(ctx.getQueryTranscriptFeaturesMax());

        transcriptComparisonResult.setTargetStart(ctx.getTargetTranscriptFeaturesMin());
        transcriptComparisonResult.setTargetStop(ctx.getTargetTranscriptFeaturesMax());

        transcriptComparisonResult.setContig(targetTranscript.getBaseData().getContig());

        transcriptComparisonResult.setQueryForwardStrand(ctx.getQueryForwardStrand());
        transcriptComparisonResult.setTargetForwardStrand(ctx.getTargetForwardStrand());
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
            if (!config.isEnabled(type)) continue;
            map.computeIfAbsent(type, k -> new ArrayList<>()).add(f);
        }
        return map;
    }

    public void compareMappedFeaturePairs(String featureType, List<FeaturePair> pairs, TranscriptComparisonResult txRes, ComparisonResult geneRes) {
        var featRes = new FeatureComparisonResult();
        featRes.setFeatureType(featureType);
        txRes.addFeatureComparison(featRes);

        if (pairs.isEmpty()) {
            featRes.setAreSameFeatures(false);
            featRes.setMissingInQueryTranscript(true);
            txRes.setAreSameTranscript(false);
            geneRes.setAreSameGene(false);
            return;
        }

        for (var p : pairs) {
            compareFeaturePair(p, featRes, txRes, geneRes);
        }
    }

    private void compareFeatures(List<FeaturePair> exonPairs, TranscriptFeature targetTranscript, TranscriptFeature queryTranscript, TranscriptComparisonResult transcriptComparisonResult, ComparisonResult geneResult, int padBp) {

        if(config.isEnabled("exon"))
            compareMappedFeaturePairs("exon", exonPairs, transcriptComparisonResult, geneResult);

        if(config.isEnabled("intron")){
            var intronPairs = mapIntronsByExonPairs(targetTranscript, queryTranscript, exonPairs);

            if(intronPairs.isEmpty()){
                boolean hadTarget = !featuresOfType(targetTranscript, "intron").isEmpty();
                boolean hadQuery = !featuresOfType(queryTranscript, "intron").isEmpty();

                if(hadTarget || hadQuery){
                    var featureComparisonResult = new FeatureComparisonResult();
                    featureComparisonResult.setFeatureType("intron");
                    featureComparisonResult.setAreSameFeatures(false);
                    if (hadTarget && !hadQuery) featureComparisonResult.setMissingInQueryTranscript(true);
                    if (!hadTarget && hadQuery) featureComparisonResult.setMissingInTargetTranscript(true);
                    transcriptComparisonResult.addFeatureComparison(featureComparisonResult);
                    transcriptComparisonResult.setAreSameTranscript(false);
                    geneResult.setAreSameGene(false);
                }
            } else{
                compareMappedFeaturePairs("intron", intronPairs, transcriptComparisonResult, geneResult);
            }
        }

        for (String ft : List.of("CDS", "UTR5", "UTR3", "start_codon", "stop_codon")) {
            if (!config.isEnabled(ft)) continue;
            var pairs = mapFeaturesWithinExonPairs(targetTranscript, queryTranscript, exonPairs, ft, padBp);
            if (pairs.isEmpty()) {
                boolean hadTarget = !featuresOfType(targetTranscript, ft).isEmpty();
                boolean hadQuery = !featuresOfType(queryTranscript, ft).isEmpty();
                if (hadTarget || hadQuery) {
                    var featureComparisonResult = new FeatureComparisonResult();
                    featureComparisonResult.setFeatureType(ft);
                    featureComparisonResult.setAreSameFeatures(false);
                    if (hadTarget && !hadQuery) featureComparisonResult.setMissingInQueryTranscript(true);
                    if (!hadTarget && hadQuery) featureComparisonResult.setMissingInTargetTranscript(true);
                    transcriptComparisonResult.addFeatureComparison(featureComparisonResult);
                    transcriptComparisonResult.setAreSameTranscript(false);
                    geneResult.setAreSameGene(false);
                }
                continue;
            }
            compareMappedFeaturePairs(ft, pairs, transcriptComparisonResult, geneResult);
        }
    }

    public static List<FeaturePair> mapIntronsByExonPairs(
            TranscriptFeature targetTranscript,
            TranscriptFeature queryTranscript,
            List<FeaturePair> exonPairs
    ) {
        var result = new ArrayList<FeaturePair>();
        var usedT = new HashSet<GtfFeature>();
        var usedQ = new HashSet<GtfFeature>();

        var tIntrons = featuresOfType(targetTranscript, "intron");
        var qIntrons = featuresOfType(queryTranscript, "intron");

        var tExonsOrdered = sortedExons(targetTranscript);
        var qExonsOrdered = sortedExons(queryTranscript);

        var tIndex = new HashMap<GtfFeature, Integer>();
        for (int i = 0; i < tExonsOrdered.size(); i++) tIndex.put(tExonsOrdered.get(i), i);

        var qIndex = new HashMap<GtfFeature, Integer>();
        for (int i = 0; i < qExonsOrdered.size(); i++) qIndex.put(qExonsOrdered.get(i), i);

        // Build direct mappings for exons (only fully mapped exon pairs).
        var t2q = new HashMap<GtfFeature, GtfFeature>();
        var q2t = new HashMap<GtfFeature, GtfFeature>();
        for (var fp : exonPairs) {
            var te = fp.getTarget();
            var qe = fp.getQuery();
            if (te != null && qe != null) {
                t2q.put(te, qe);
                q2t.put(qe, te);
            }
        }

        for (int i = 0; i + 1 < tExonsOrdered.size(); i++) {
            var teL = tExonsOrdered.get(i);
            var teR = tExonsOrdered.get(i + 1);

            var qeL = t2q.get(teL);
            var qeR = t2q.get(teR);
            if (qeL == null || qeR == null) {
                continue;
            }

            var qiL = qIndex.get(qeL);
            var qiR = qIndex.get(qeR);
            if (qiL == null || qiR == null) continue;

            // Require adjacency and correct order on the query side as well (no crossing).
            if (qiR != qiL + 1) {
                // Not strictly adjacent on query; skip pairing for this intron to avoid ambiguity.
                continue;
            }

            // Find the introns bounded by the exon pairs on both sides.
            var tIntron = findIntronBetween(teL, teR, tIntrons);
            var qIntron = findIntronBetween(qeL, qeR, qIntrons);

            if (tIntron != null || qIntron != null) {
                result.add(new FeaturePair(tIntron, qIntron));
                if (tIntron != null) usedT.add(tIntron);
                if (qIntron != null) usedQ.add(qIntron);
            }
        }

        for (var ti : tIntrons) {
            if (!usedT.contains(ti)) result.add(new FeaturePair(ti, null));
        }
        for (var qi : qIntrons) {
            if (!usedQ.contains(qi)) result.add(new FeaturePair(null, qi));
        }

        return result;
    }

    private static GtfFeature findIntronBetween(
            GtfFeature eLeft,
            GtfFeature eRight,
            List<GtfFeature> introns
    ) {
        int exonLeftEnd = eLeft.getBaseData().getEnd() + 1;
        int exonRightStart = eRight.getBaseData().getStart() - 1;
        if (exonRightStart < exonLeftEnd) return null;

        GtfFeature best = null;
        int bestOverlap = -1;

        for (var intron : introns) {
            int start = intron.getBaseData().getStart();
            int end = intron.getBaseData().getEnd();

            if (!within(start, end, exonLeftEnd, exonRightStart, 0)) continue;

            int overlap = Math.min(end, exonRightStart) - Math.max(start, exonLeftEnd) + 1;
            if (overlap > bestOverlap) {
                bestOverlap = overlap;
                best = intron;
            }
        }
        return best;
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

        var targetTranscriptBiotype = transcriptComparisonResult.getTargetBiotype();
        var queryTranscriptBiotype = transcriptComparisonResult.getQueryBiotype();

        var ctx = new ComparisonContext(targetFeature, queryFeature, targetTranscriptBiotype, queryTranscriptBiotype, config, targetSequenceExtractor, querySequenceExtractor, null, null);
        var currentFeatureType = targetBaseData.getType();
        ServiceLoader<? extends ComparisonFeature> currentLoader = featureLoader;
        if (Constants.CDS.equals(GtfConfig.getDefault(currentFeatureType))) {
            currentLoader = cdsLoader;
        }
        for (var comp : currentLoader) {
            if (!config.isEnabled(comp.getName())) continue;

            boolean changed;
            try {
                changed = comp.compare(ctx);
            } catch (Exception e) {
                throw new RuntimeException(e);
            }

            if (changed) {
                addToRegionComparison(comp.getName(), regionComparisonResult, featureComparisonResult, transcriptComparisonResult, geneResult);
            }

            logger.debug("Feature {}: {} changed: {}", targetFeature.getBaseData().getType(), comp.getName(), changed);
        }
    }

    private void addToRegionComparison(String name, RegionComparisonResult regionComparisonResult, FeatureComparisonResult featureComparisonResult, TranscriptComparisonResult transcriptComparisonResult, ComparisonResult result) {
        var geneComparisonResult = result.getGeneComparison();
        switch (name) {
            case Constants.LENGTH_COMPARATOR_NAME:
                transcriptComparisonResult.setAreSameTranscript(false);
                result.setAreSameGene(false);
                geneComparisonResult.setAreSameGene(false);
                featureComparisonResult.setAreSameFeatures(false);

                regionComparisonResult.setLengthDifferent(true);
                break;
            case Constants.SEQUENCE_COMPARATOR_NAME:
                transcriptComparisonResult.setAreSameTranscript(false);
                geneComparisonResult.setAreSameGene(false);
                result.setAreSameGene(false);
                featureComparisonResult.setAreSameFeatures(false);

                regionComparisonResult.setSequenceDifferenceFound(true);
                transcriptComparisonResult.setSequenceDifferent(true);
                geneComparisonResult.setSequenceDifferent(true);
                // TODO Add Sequence Difference
                break;
            case Constants.START_COMPARATOR_NAME:
                transcriptComparisonResult.setAreSameTranscript(false);
                geneComparisonResult.setAreSameGene(false);
                result.setAreSameGene(false);
                featureComparisonResult.setAreSameFeatures(false);

                regionComparisonResult.setStartDifferent(true);
                break;
            case Constants.STOP_COMPARATOR_NAME:
                transcriptComparisonResult.setAreSameTranscript(false);
                geneComparisonResult.setAreSameGene(false);
                result.setAreSameGene(false);
                featureComparisonResult.setAreSameFeatures(false);

                regionComparisonResult.setEndDifferent(true);
                break;
            case Constants.SAME_PROTEIN_COMPARATOR_NAME:
                transcriptComparisonResult.setAreSameTranscript(false);
                geneComparisonResult.setAreSameGene(false);
                result.setAreSameGene(false);
                featureComparisonResult.setAreSameFeatures(false);

                regionComparisonResult.setProteinDifferent(true);
                break;
            default:
                logger.warn("Unknown gene comparison feature: {}", name);
                break;
        }
    }

    private void addToTranscriptComparison(String name, TranscriptComparisonResult transcriptComparisonResult, ComparisonResult result) {
        var geneComparisonResult = result.getGeneComparison();
        switch (name) {
            case Constants.TRANSCRIPT_LENGTH_COMPARATOR_NAME:
                transcriptComparisonResult.setAreSameTranscript(false);
                result.setAreSameGene(false);
                geneComparisonResult.setAreSameGene(false);

                transcriptComparisonResult.setLengthDifferent(true);
                break;
            case Constants.TRANSCRIPT_START_COMPARATOR_NAME:
                transcriptComparisonResult.setAreSameTranscript(false);
                geneComparisonResult.setAreSameGene(false);
                result.setAreSameGene(false);

                transcriptComparisonResult.setStartDifferent(true);

                break;
            case Constants.TRANSCRIPT_STOP_COMPARATOR_NAME:
                transcriptComparisonResult.setAreSameTranscript(false);
                geneComparisonResult.setAreSameGene(false);
                result.setAreSameGene(false);

                transcriptComparisonResult.setStopDifferent(true);
                break;
            case Constants.TRANSCRIPT_SEQUENCE_COMPARATOR_NAME:
                transcriptComparisonResult.setAreSameTranscript(false);
                geneComparisonResult.setAreSameGene(false);
                result.setAreSameGene(false);

                transcriptComparisonResult.setStrandDifferent(true);
                break;

            case Constants.TRANSCRIPT_BIOTYPE_COMPARATOR_NAME:
                transcriptComparisonResult.setAreSameTranscript(false);
                geneComparisonResult.setAreSameGene(false);
                result.setAreSameGene(false);

                transcriptComparisonResult.setBiotypeDifferent(true);
                break;
            case Constants.TRANSCRIPT_STRAND_COMPARATOR_NAME:
                transcriptComparisonResult.setAreSameTranscript(false);
                geneComparisonResult.setAreSameGene(false);
                result.setAreSameGene(false);

                transcriptComparisonResult.setStrandDifferent(true);
                break;
            default:
                logger.warn("Unknown transcript comparison feature: {}", name);
                break;
        }
    }

    private void addToGeneComparison(String name, ComparisonResult result) {
        var geneComparisonResult = result.getGeneComparison();
        switch (name) {
            case Constants.GENE_LENGTH_COMPARATOR_NAME:
                result.setAreSameGene(false);
                geneComparisonResult.setAreSameGene(false);

                geneComparisonResult.setLengthDifferent(true);
                break;
            case Constants.GENE_START_COMPARATOR_NAME:
                geneComparisonResult.setAreSameGene(false);
                result.setAreSameGene(false);

                geneComparisonResult.setStartDifferent(true);

                break;
            case Constants.GENE_STOP_COMPARATOR_NAME:
                geneComparisonResult.setAreSameGene(false);
                result.setAreSameGene(false);

                geneComparisonResult.setStopDifferent(true);
                break;
            case Constants.GENE_STRAND_COMPARATOR_NAME:
                geneComparisonResult.setAreSameGene(false);
                result.setAreSameGene(false);

                geneComparisonResult.setStrandDifferent(true);
                break;

            case Constants.GENE_CONTIG_COMPARATOR_NAME:
                geneComparisonResult.setAreSameGene(false);
                result.setAreSameGene(false);

                geneComparisonResult.setContigDifferent(true);
                break;
            default:
                logger.warn("Unknown gene comparison feature: {}", name);
                break;
        }
    }

    private void markMissingFeature(GtfFeature targetFeature, GtfFeature queryFeature, FeatureComparisonResult featRes, TranscriptComparisonResult transcriptComparisonResult, ComparisonResult comparisonResult) {

        comparisonResult.setAreSameGene(false);
        comparisonResult.getGeneComparison().setAreSameGene(false);
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
