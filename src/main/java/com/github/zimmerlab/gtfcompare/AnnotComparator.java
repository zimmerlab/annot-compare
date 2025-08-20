package com.github.zimmerlab.gtfcompare;


import com.github.kleinsamuel.gtfutils.GtfConfig;
import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.GeneFeature;
import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.compare.*;
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

        ResultWriter.writeComparisonResult(comparisonResults, outputPath);
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

        var pairedExons = pairExonsByGapAlignment(targetTranscript, queryTranscript, 2, 200, 50, 0.1);

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

    public void compareMappedFeaturePairs(
            String featureType,
            List<FeaturePair> pairs,
            TranscriptComparisonResult txRes,
            ComparisonResult geneRes
    ) {
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
            compareFeaturePair(p, featRes, txRes, geneRes); // deine bestehende Methode
        }
    }

    private void compareFeatures(List<FeaturePair> exonPairs, TranscriptFeature targetTranscript, TranscriptFeature queryTranscript, TranscriptComparisonResult transcriptComparisonResult, ComparisonResult geneResult, int padBp) {

        compareMappedFeaturePairs("exon", exonPairs, transcriptComparisonResult, geneResult);

        for (String ft : List.of("CDS","UTR5","UTR3","start_codon","stop_codon")) {
            var pairs = mapFeaturesWithinExonPairs(targetTranscript, queryTranscript, exonPairs, ft, padBp);
            if (pairs.isEmpty()) {
                // prüfen, ob es auf einer Seite überhaupt Features dieses Typs gab
                boolean hadTarget = !featuresOfType(targetTranscript, ft).isEmpty();
                boolean hadQuery  = !featuresOfType(queryTranscript, ft).isEmpty();
                if (hadTarget || hadQuery) {
                    var fr = new FeatureComparisonResult();
                    fr.setFeatureType(ft);
                    fr.setAreSameFeatures(false);
                    if (hadTarget && !hadQuery) fr.setMissingInQueryTranscript(true);
                    if (!hadTarget && hadQuery) fr.setMissingInTargetTranscript(true);
                    transcriptComparisonResult.addFeatureComparison(fr);
                    transcriptComparisonResult.setAreSameTranscript(false);
                    geneResult.setAreSameGene(false);
                }
                continue;
            }
            compareMappedFeaturePairs(ft, pairs, transcriptComparisonResult, geneResult);
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


    private List<FeaturePair> pairByPosition(List<GtfFeature> targets, List<GtfFeature> queries) {

        Comparator<GtfFeature> byStartThenEnd = Comparator.comparingInt((GtfFeature f) -> f.getBaseData().getStart()).thenComparingInt(f -> f.getBaseData().getEnd());

        var a = new ArrayList<>(targets);
        var b = new ArrayList<>(queries);
        a.sort(byStartThenEnd);
        b.sort(byStartThenEnd);

        var aSize = a.size();
        var bSize = b.size();
        var n = Math.max(aSize, bSize);
        var pairs = new ArrayList<FeaturePair>(n);
        for (int i = 0; i < n; i++) {
            GtfFeature ta = i < aSize ? a.get(i) : null;
            GtfFeature qb = i < bSize ? b.get(i) : null;
            pairs.add(new FeaturePair(ta, qb));
        }
        return pairs;
    }

    private List<FeaturePair> pairByOverlap(List<GtfFeature> targets, List<GtfFeature> queries) {
        var pairs = new ArrayList<FeaturePair>();
        for (GtfFeature ta : targets) {
            if (ta.getBaseData().getAttributes("transcript_id") != null && ta.getBaseData().getAttributes("transcript_id").getFirst().equals("ENST00000412513")) {
                var a = 2;
            }
            var targetStart = ta.getBaseData().getStart();
            var targetEnd = ta.getBaseData().getEnd();

            var found = false;
            for (GtfFeature qb : queries) {
                var queryStart = qb.getBaseData().getStart();
                int queryEnd = qb.getBaseData().getEnd();

                if (queryStart <= targetEnd && queryEnd >= targetStart) {
                    pairs.add(new FeaturePair(ta, qb));
                    found = true;
                }
            }

            if (!found) {
                pairs.add(new FeaturePair(ta, null));
            }
        }

        // Option: Queries, die nie gematcht wurden, mit null target zurückgeben
        for (GtfFeature qb : queries) {
            var matched = pairs.stream().anyMatch(p -> p.getQuery() != null && p.getQuery().equals(qb));
            if (!matched) {
                pairs.add(new FeaturePair(null, qb));
            }
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


    private static List<GtfFeature> sortedExons(TranscriptFeature tf) {
        boolean fwd = tf.getBaseData().isForwardStrand();
        Comparator<GtfFeature> byStart = Comparator.comparingInt(f -> f.getBaseData().getStart());
        var list = tf.getFeatures().stream().filter(f -> "exon".equals(f.getBaseData().getType())).sorted(byStart).toList();
        return fwd ? list : new ArrayList<>(list) {{
            Collections.reverse(this);
        }};
    }

    private static List<Integer> gapProfile(List<GtfFeature> exons) {
        var gaps = new ArrayList<Integer>(Math.max(0, exons.size() - 1));
        for (int i = 0; i < exons.size() - 1; i++) {
            int g = exons.get(i + 1).getBaseData().getStart() - exons.get(i).getBaseData().getEnd() - 1;
            gaps.add(g);
        }
        return gaps;
    }

    private static int[][] nwBacktrace(int[] A, int[] B, int gapPenalty, int cap) {
        int n = A.length, m = B.length;
        int[][] dp = new int[n + 1][m + 1];
        int[][] bt = new int[n + 1][m + 1]; // 0=diag,1=up,2=left

        for (int i = 1; i <= n; i++) {
            dp[i][0] = dp[i - 1][0] - gapPenalty;
            bt[i][0] = 1;
        }
        for (int j = 1; j <= m; j++) {
            dp[0][j] = dp[0][j - 1] - gapPenalty;
            bt[0][j] = 2;
        }

        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                int diff = Math.abs(A[i - 1] - B[j - 1]);
                int match = dp[i - 1][j - 1] - Math.min(diff, cap);
                int delA = dp[i - 1][j] - gapPenalty;
                int delB = dp[i][j - 1] - gapPenalty;

                if (match >= delA && match >= delB) {
                    bt[i][j] = 0;
                    dp[i][j] = match;
                } else if (delA >= delB) {
                    bt[i][j] = 1;
                    dp[i][j] = delA;
                } else {
                    bt[i][j] = 2;
                    dp[i][j] = delB;
                }
            }
        }
        return bt;
    }

    private static List<int[]> recoverGapMatches(int[][] bt, int n, int m) {
        var pairs = new ArrayList<int[]>();
        int i = n, j = m;
        while (i > 0 || j > 0) {
            if (i > 0 && j > 0 && bt[i][j] == 0) {
                pairs.add(new int[]{i - 1, j - 1});
                i--;
                j--;
            } else if (i > 0 && (j == 0 || bt[i][j] == 1)) {
                i--;
            } else {
                j--;
            }
        }
        Collections.reverse(pairs);
        return pairs;
    }

    private static int endpointManhattan(GtfFeature a, GtfFeature b) {
        return Math.abs(a.getBaseData().getStart() - b.getBaseData().getStart()) + Math.abs(a.getBaseData().getEnd() - b.getBaseData().getEnd());
    }

    public static List<FeaturePair> pairExonsByGapAlignment(
            TranscriptFeature t, TranscriptFeature q,
            int gapPenalty, int capDelta,
            int lenCapBp, double lenCapFrac
    ) {
        var pairs = new ArrayList<FeaturePair>();
        var targetExons = sortedExons(t);
        var queryExons  = sortedExons(q);

        if (targetExons.isEmpty() && queryExons.isEmpty()) return pairs;
        if (targetExons.size() == 1 && queryExons.size() == 1) {
            pairs.add(new FeaturePair(targetExons.getFirst(), queryExons.getFirst()));
            return pairs;
        }

        var A = gapProfile(targetExons).stream().mapToInt(x -> x).toArray();
        var B = gapProfile(queryExons).stream().mapToInt(x -> x).toArray();

        java.util.function.ToIntFunction<GtfFeature> len = f ->
                f.getBaseData().getEnd() - f.getBaseData().getStart() + 1;

        java.util.function.BiPredicate<GtfFeature,GtfFeature> passLen =
                (te, qe) -> {
                    int lt = len.applyAsInt(te), lq = len.applyAsInt(qe);
                    int abs = Math.abs(lt - lq);
                    double rel = (double) abs / Math.max(lt, lq);
                    return abs <= lenCapBp || rel <= lenCapFrac;
                };

        if (A.length == 0 || B.length == 0) {
            var usedQ = new HashSet<GtfFeature>();
            for (var te : targetExons) {
                GtfFeature best = null;
                int bestLenDiff = Integer.MAX_VALUE;
                int bestD = Integer.MAX_VALUE;
                for (var qe : queryExons) {
                    if (usedQ.contains(qe)) continue;
                    if (!passLen.test(te, qe)) continue;
                    int lenDiff = Math.abs(len.applyAsInt(te) - len.applyAsInt(qe));
                    int d = endpointManhattan(te, qe);
                    if (lenDiff < bestLenDiff || (lenDiff == bestLenDiff && d < bestD)) {
                        bestLenDiff = lenDiff; bestD = d; best = qe;
                    }
                }
                pairs.add(new FeaturePair(te, best));
                if (best != null) usedQ.add(best);
            }
            for (var qe : queryExons)
                if (pairs.stream().noneMatch(p -> qe.equals(p.getQuery())))
                    pairs.add(new FeaturePair(null, qe));
            return pairs;
        }

        var bt = nwBacktrace(A, B, gapPenalty, capDelta);
        var gapMatches = recoverGapMatches(bt, A.length, B.length);

        var usedT = new boolean[targetExons.size()];
        var usedQ = new boolean[queryExons.size()];

        var deltas = new ArrayList<Integer>();
        java.util.function.IntSupplier deltaHat = () -> {
            if (deltas.size() < 2) return 0;
            var copy = new ArrayList<>(deltas);
            copy.sort(Integer::compare);
            int mid = copy.size()/2;
            return (copy.size()%2==1) ? copy.get(mid) : (copy.get(mid-1)+copy.get(mid))/2;
        };

        for (var g : gapMatches) {
            int k = g[0], l = g[1];
            int[] candT = {k, k + 1};
            int[] candQ = {l, l + 1};

            record Cand(int ct, int cq, int lenDiff, int rank) {}
            var cands = new ArrayList<Cand>();

            int Δ = deltaHat.getAsInt();

            for (int ct : candT)
                for (int cq : candQ) {
                    if (ct < 0 || cq < 0 || ct >= targetExons.size() || cq >= queryExons.size()) continue;
                    if (usedT[ct] || usedQ[cq]) continue;
                    var te = targetExons.get(ct);
                    var qe = queryExons.get(cq);
                    if (!passLen.test(te, qe)) continue;

                    int lenDiff = Math.abs(len.applyAsInt(te) - len.applyAsInt(qe));

                    int rank = Math.abs((te.getBaseData().getStart() - qe.getBaseData().getStart()) - Δ)
                            + Math.abs((te.getBaseData().getEnd()   - qe.getBaseData().getEnd())   - Δ);

                    cands.add(new Cand(ct, cq, lenDiff, rank));
                }

            cands.sort((x,y) -> {
                if (x.lenDiff != y.lenDiff) return Integer.compare(x.lenDiff, y.lenDiff); // 1) min lenDiff
                return Integer.compare(x.rank, y.rank);                                    // 2) min rank
            });

            if (!cands.isEmpty()) {
                var best = cands.get(0);
                var te = targetExons.get(best.ct);
                var qe = queryExons.get(best.cq);
                pairs.add(new FeaturePair(te, qe));
                usedT[best.ct] = usedQ[best.cq] = true;

                deltas.add(te.getBaseData().getStart() - qe.getBaseData().getStart());
            }
        }

        // left anchor
        int iL = 0, jL = 0;
        while (iL < targetExons.size() && usedT[iL]) iL++;
        while (jL < queryExons.size()  && usedQ[jL]) jL++;
        if (iL < targetExons.size() && jL < queryExons.size()) {
            var te = targetExons.get(iL);
            var qe = queryExons.get(jL);
            int lt = te.getBaseData().getEnd() - te.getBaseData().getStart() + 1;
            int lq = qe.getBaseData().getEnd() - qe.getBaseData().getStart() + 1;
            int lenDiff = Math.abs(lt - lq);
            double rel = (double) lenDiff / Math.max(lt, lq);
            if (lenDiff <= lenCapBp || rel <= lenCapFrac) {
                pairs.add(new FeaturePair(te, qe));
                usedT[iL] = usedQ[jL] = true;
            }
        }

// right anchor
        int iR = targetExons.size() - 1, jR = queryExons.size() - 1;
        while (iR >= 0 && usedT[iR]) iR--;
        while (jR >= 0 && usedQ[jR]) jR--;
        if (iR >= 0 && jR >= 0) {
            var te = targetExons.get(iR);
            var qe = queryExons.get(jR);
            int lt = te.getBaseData().getEnd() - te.getBaseData().getStart() + 1;
            int lq = qe.getBaseData().getEnd() - qe.getBaseData().getStart() + 1;
            int lenDiff = Math.abs(lt - lq);
            double rel = (double) lenDiff / Math.max(lt, lq);
            if (lenDiff <= lenCapBp || rel <= lenCapFrac) {
                pairs.add(new FeaturePair(te, qe));
                usedT[iR] = usedQ[jR] = true;
            }
        }

// greedy for the rest
        int Δ = 0;
        if (!pairs.isEmpty()) {
            var ds = new ArrayList<Integer>();
            for (var p : pairs) {
                if (p.getTarget()!=null && p.getQuery()!=null) {
                    ds.add(p.getTarget().getBaseData().getStart() - p.getQuery().getBaseData().getStart());
                }
            }
            if (ds.size() >= 2) {
                ds.sort(Integer::compare);
                int mid = ds.size()/2;
                Δ = (ds.size()%2==1) ? ds.get(mid) : (ds.get(mid-1)+ds.get(mid))/2;
            }
        }
        for (int it = 0; it < targetExons.size(); it++) {
            if (usedT[it]) continue;
            var te = targetExons.get(it);
            GtfFeature best = null;
            int bestLenDiff = Integer.MAX_VALUE;
            int bestRank   = Integer.MAX_VALUE;
            for (int iq = 0; iq < queryExons.size(); iq++) {
                if (usedQ[iq]) continue;
                var qe = queryExons.get(iq);
                int lt = te.getBaseData().getEnd() - te.getBaseData().getStart() + 1;
                int lq = qe.getBaseData().getEnd() - qe.getBaseData().getStart() + 1;
                int lenDiff = Math.abs(lt - lq);
                double rel  = (double) lenDiff / Math.max(lt, lq);
                if (!(lenDiff <= lenCapBp || rel <= lenCapFrac)) continue;

                int rank = Math.abs((te.getBaseData().getStart() - qe.getBaseData().getStart()) - Δ)
                        + Math.abs((te.getBaseData().getEnd()   - qe.getBaseData().getEnd())   - Δ);

                if (lenDiff < bestLenDiff || (lenDiff == bestLenDiff && rank < bestRank)) {
                    bestLenDiff = lenDiff; bestRank = rank; best = qe;
                }
            }
            pairs.add(new FeaturePair(te, best));
            if (best != null) {
                int idx = queryExons.indexOf(best);
                usedT[it] = true; usedQ[idx] = true;
            }
        }

        for (int i = 0; i < targetExons.size(); i++)
            if (!usedT[i]) pairs.add(new FeaturePair(targetExons.get(i), null));
        for (int j = 0; j < queryExons.size(); j++)
            if (!usedQ[j]) pairs.add(new FeaturePair(null, queryExons.get(j)));

        return pairs;
    }


    private static boolean within(int s, int e, int xs, int xe, int pad) {
        return s >= xs - pad && e <= xe + pad;
    }

    private static List<GtfFeature> featuresOfType(TranscriptFeature tf, String type) {
        return tf.getFeatures().stream()
                .filter(f -> type.equals(com.github.kleinsamuel.gtfutils.GtfConfig.getDefault(f.getBaseData().getType())))
                .sorted(Comparator.comparingInt(f -> f.getBaseData().getStart()))
                .toList();
    }

    public static List<FeaturePair> mapFeaturesWithinExonPairs(
            TranscriptFeature t,
            TranscriptFeature q,
            List<FeaturePair> exonPairs,      // aus pairExonsByGapAlignment(t,q,...)
            String featureType,
            int padBp                          // z.B. 10..50 bp
    ) {
        var res = new ArrayList<FeaturePair>();
        var usedQ = new HashSet<GtfFeature>();

        var tFeat = featuresOfType(t, featureType);
        var qFeat = featuresOfType(q, featureType);

        // Schneller Zugriff: Query-Features pro Query-Exon sammeln
        var qByExon = new HashMap<GtfFeature, List<GtfFeature>>();
        for (var fp : exonPairs) {
            var qe = fp.getQuery();
            if (qe != null) qByExon.put(qe, new ArrayList<>());
        }
        for (var qf : qFeat) {
            // finde Query-Exon-Container (erstes passendes Exon reicht)
            for (var fp : exonPairs) {
                var qe = fp.getQuery();
                if (qe == null) continue;
                var qS = qf.getBaseData().getStart();
                var qE = qf.getBaseData().getEnd();
                var qeS = fp.getQuery().getBaseData().getStart();
                var qeE = fp.getQuery().getBaseData().getEnd();
                if (within(qS, qE, qeS, qeE, padBp)) {
                    qByExon.computeIfAbsent(qe, k -> new ArrayList<>()).add(qf);
                    break;
                }
            }
        }

        for (var fpExon : exonPairs) {
            var te = fpExon.getTarget();
            var qe = fpExon.getQuery();
            // alle Target-Features, die in dieses Target-Exon fallen
            var tSub = new ArrayList<GtfFeature>();
            if (te != null) {
                int teS = te.getBaseData().getStart(), teE = te.getBaseData().getEnd();
                for (var tf : tFeat) {
                    int s = tf.getBaseData().getStart(), e = tf.getBaseData().getEnd();
                    if (within(s, e, teS, teE, padBp)) tSub.add(tf);
                }
            }

            var candidates = (qe != null) ? qByExon.getOrDefault(qe, List.of()) : List.<GtfFeature>of();

            for (var tf : tSub) {
                GtfFeature best = null;
                int bestDist = Integer.MAX_VALUE;
                for (var qf : candidates) {
                    if (usedQ.contains(qf)) continue;
                    int d = Math.abs(tf.getBaseData().getStart() - qf.getBaseData().getStart())
                            + Math.abs(tf.getBaseData().getEnd()   - qf.getBaseData().getEnd());
                    if (d < bestDist) { bestDist = d; best = qf; }
                }
                res.add(new FeaturePair(tf, best));
                if (best != null) usedQ.add(best);
            }
        }

        for (var qf : qFeat) {
            if (!usedQ.contains(qf) && res.stream().noneMatch(p -> qf.equals(p.getQuery()))) {
                res.add(new FeaturePair(null, qf));
            }
        }
        return res;
    }
}
