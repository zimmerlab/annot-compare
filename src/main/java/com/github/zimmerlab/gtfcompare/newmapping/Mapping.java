package com.github.zimmerlab.gtfcompare.newmapping;

import com.github.kleinsamuel.gtfutils.GtfConfig;
import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.model.GenePair;
import com.github.zimmerlab.gtfcompare.utils.Constants;

import java.util.*;
import java.util.stream.Collectors;

public class Mapping {
    private record CigarOp(String type, int length) {
    }

    private record LeftoverGeneIds(Set<String> targetOnly, Set<String> queryOnly) {
    }

    private static final double SIMILARITY_CUTOFF = 0.99;
    private static final int SHORT_EXON_THRESHOLD = 3000;
    private static final double STRICT_SHORT_MIN_SIM = 1.0;

    public static List<GenePair> map(GtfFile targetGtf, GtfFile queryGtf, Set<String> allowedTypes) {
        var targetGeneIds = new HashSet<>(targetGtf.getAllGeneFeatureIds());
        var queryGeneIds = new HashSet<>(queryGtf.getAllGeneFeatureIds());

        var idMapping = idMapping(targetGtf, queryGtf, targetGeneIds, queryGeneIds);
        var noSameTranscripts = getGenesNoSameTranscript(idMapping);
        var leftOvers = getLeftoverGeneIds(targetGeneIds, queryGeneIds);

        if (allowedTypes == null) allowedTypes = new HashSet<>(Constants.FEATURE_TYPES);

        Set<String> finalAllowedTypes = allowedTypes;
        var targetLeftOverCigars = leftOvers.targetOnly
                .stream()
                .collect(Collectors.toMap(
                        geneId -> geneId,
                        geneId -> targetGtf.getGeneFeature(geneId)
                                .getTranscripts()
                                .stream()
                                .map(transcript -> structureCigar(transcript, finalAllowedTypes))
                                .collect(Collectors.toList())
                ));

        var queryLeftOverCigars = leftOvers.queryOnly
                .stream()
                .collect(Collectors.toMap(
                        geneId -> geneId,
                        geneId -> queryGtf.getGeneFeature(geneId)
                                .getTranscripts()
                                .stream()
                                .map(transcript -> structureCigar(transcript, finalAllowedTypes))
                                .collect(Collectors.toList())
                ));

        var results = new LinkedList<GenePair>();
        var allQueryGenes  = new HashSet<>(queryLeftOverCigars.keySet());
        var allTargetGenes = new HashSet<>(targetLeftOverCigars.keySet());

        var noSameTranscriptsLeftovers = new HashSet<GenePair>();
        for (var genePair : noSameTranscripts) {
            var target = genePair.getTarget();
            var query = genePair.getQuery();

            var bestScore = Double.NEGATIVE_INFINITY;

            for(var targetTranscript : target.getTranscripts()){
                var targetCigar = structureCigar(targetTranscript, allowedTypes);

                for(var queryTranscript : query.getTranscripts()){
                    var queryCigar = structureCigar(queryTranscript, allowedTypes);

                    var score = cigarSimilarity(targetCigar, queryCigar);

                    if(score > bestScore){
                        bestScore = score;
                    }
                }
            }

            if(bestScore > SIMILARITY_CUTOFF){
                results.add(genePair);
            } else{
                noSameTranscriptsLeftovers.add(genePair);
            }

        }

        for (var leftOverQueryCigars : queryLeftOverCigars.entrySet()) {
            var queryGeneId = leftOverQueryCigars.getKey();
            var queryCigars = leftOverQueryCigars.getValue();
            for (var leftOverTargetCigars : targetLeftOverCigars.entrySet()) {
                var targetGeneId = leftOverTargetCigars.getKey();
                var targetCigars = leftOverTargetCigars.getValue();
                var bestScore = Double.NEGATIVE_INFINITY;
                if(targetGeneId.equals("ENSG00000274265") && queryGeneId.equals("ENSG00000269501")){
                    var a = 2;
                }
                for (var targetCigar : targetCigars) {
                    for (var queryCigar : queryCigars) {
                        var score = cigarSimilarity(targetCigar, queryCigar);
                        bestScore = Math.max(bestScore, score);

                        if(score > SIMILARITY_CUTOFF){
                            var a = 2;
                        }
                    }
                }
                if (bestScore > SIMILARITY_CUTOFF) {
                    results.add(new GenePair(targetGtf.getGeneFeature(targetGeneId), queryGtf.getGeneFeature(queryGeneId)));

                    allQueryGenes.remove(queryGeneId);
                    allTargetGenes.remove(targetGeneId);
                }

            }
        }

        return null;
    }

    private static List<GenePair> idMapping(GtfFile targetGtf, GtfFile queryGtf, Set<String> targetGeneIds, Set<String> queryGeneIds) {
        var commonIds = new HashSet<>(targetGeneIds);
        commonIds.retainAll(queryGeneIds);
        var mapping = new ArrayList<GenePair>(commonIds.size());

        for (var geneId : commonIds) {
            mapping.add(new GenePair(targetGtf.getGeneFeature(geneId), queryGtf.getGeneFeature(geneId)));
        }
        return mapping;
    }

    private static List<GenePair> getGenesNoSameTranscript(List<GenePair> idMapping) {
        var noSameTranscripts = new ArrayList<GenePair>();
        for (var genePair : idMapping) {
            var target = genePair.getTarget();
            var query = genePair.getQuery();
            var hasSameTranscript = false;

            outerLoop:
            for (var targetTranscript : target.getTranscripts()) {
                var targetTranscriptId = targetTranscript.getTranscriptId();
                for (var queryTranscript : query.getTranscripts()) {
                    var queryTranscriptId = queryTranscript.getTranscriptId();

                    if (queryTranscriptId.equals(targetTranscriptId)) {
                        hasSameTranscript = true;
                        break outerLoop;
                    }
                }
            }

            if (!hasSameTranscript) noSameTranscripts.add(genePair);
        }
        return noSameTranscripts;
    }


    private static LeftoverGeneIds getLeftoverGeneIds(Set<String> targetGeneIds, Set<String> queryGeneIds) {
        var targetOnly = new HashSet<>(targetGeneIds);
        targetOnly.removeAll(queryGeneIds);

        var queryOnly = new HashSet<>(queryGeneIds);
        queryOnly.removeAll(targetGeneIds);

        return new LeftoverGeneIds(targetOnly, queryOnly);
    }

    private static List<CigarOp> structureCigar(TranscriptFeature feature, Set<String> allowedTypes) {
        var result = new ArrayList<CigarOp>();

        feature.getFeatures()
                .stream()
                .sorted(Comparator.comparingInt(f -> f.getBaseData().getStart()))
                .forEach(g -> {
                    var baseData = g.getBaseData();
                    var type = baseData.getType();
                    var defaultType = GtfConfig.getDefault(type);

                    if (!allowedTypes.contains(defaultType)) return;

                    var start = baseData.getStart();
                    var stop = baseData.getEnd();

                    var len = stop - start + 1;
                    result.add(new CigarOp(type, len));
                });

        return result;
    }

    private static double cigarSimilarity(List<CigarOp> a, List<CigarOp> b) {
        if (a.isEmpty() && b.isEmpty()) return 1.0;

        if (a.isEmpty() || b.isEmpty()) return 0.0;

        if(a.size() == 1 && b.size() == 1){
            return handleSingleExon(a, b);
        }

        var minSize = Math.min(a.size(), b.size());
        var maxSize = Math.max(a.size(), b.size());

        double similarityScore = 0.0;

        for (int i = 0; i < minSize; i++) {
            var opA = a.get(i);
            var opB = b.get(i);

            if (!opA.type().equals(opB.type())) continue;

            int lenA = opA.length();
            int lenB = opB.length();

            int maxLen = Math.max(lenA, lenB);

            if (maxLen == 0) {
                similarityScore += 1;
            } else {
                var relDiff = Math.abs(lenA - lenB) / (double) maxLen;
                double posSimilarity = 1.0 - relDiff;

                if(maxLen <= SHORT_EXON_THRESHOLD){
                    posSimilarity = posSimilarity >= STRICT_SHORT_MIN_SIM ? posSimilarity : 0.0;
                }

                similarityScore += posSimilarity;
            }
        }

        return similarityScore / (double) maxSize;
    }

    private static double handleSingleExon(List<CigarOp> a, List<CigarOp> b) {
        var opA = a.getFirst();
        var opB = b.getFirst();

        if(!opA.type().equals(opB.type())) return 0.0;

        var maxLength = Math.max(opA.length(), opB.length());

        if(maxLength == 0) return 0.0;

        var relDiff = Math.abs(opA.length() - opB.length()) / (double) maxLength;
        var baseSimilarity = 1.0 - relDiff;

        if(maxLength < SHORT_EXON_THRESHOLD){
            return baseSimilarity > STRICT_SHORT_MIN_SIM ? baseSimilarity : 0.0;
        }

        return baseSimilarity / (double) maxLength;
    }
}