package com.github.zimmerlab.gtfcompare.newmapping;

import com.github.kleinsamuel.gtfutils.GtfConfig;
import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.model.GenePair;
import com.github.zimmerlab.gtfcompare.utils.Constants;

import java.util.*;
import java.util.stream.Collectors;

import static com.github.zimmerlab.gtfcompare.newmapping.Cigar.cigarSimilarity;
import static com.github.zimmerlab.gtfcompare.newmapping.Cigar.structureCigar;
import static com.github.zimmerlab.gtfcompare.newmapping.MappingConstants.SIMILARITY_CUTOFF;

// TODO single exons sequenz homology?

public class Mapping {

    private record LeftoverGeneIds(Set<String> targetOnly, Set<String> queryOnly) {
    }



    public static List<List<GenePair>> map(GtfFile targetGtf, GtfFile queryGtf, Set<String> allowedTypes) {
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
        var unmappedQueryGenes = new HashSet<>(queryLeftOverCigars.keySet());
        var unmappedTargetGenes = new HashSet<>(targetLeftOverCigars.keySet());

        var noSameTranscriptsLeftovers = new HashSet<GenePair>();

        var targetCigars = new HashMap<String, List<List<CigarOp>>>();
        var queryCigars = new HashMap<String, List<List<CigarOp>>>();
        for (var genePair : idMapping) {
            var target = genePair.getTarget();
            var query = genePair.getQuery();

            var targetGeneId = target.getGeneId();
            var queryGeneId = query.getGeneId();

            var bestScore = Double.NEGATIVE_INFINITY;

            var firstRun = true;
            for (var targetTranscript : target.getTranscripts()) {
                var targetCigar = structureCigar(targetTranscript, allowedTypes);
                targetCigars.computeIfAbsent(targetGeneId, geneId -> new LinkedList<>()).add(targetCigar);
                for (var queryTranscript : query.getTranscripts()) {
                    var queryCigar = structureCigar(queryTranscript, allowedTypes);
                    if(firstRun) {
                        queryCigars.computeIfAbsent(queryGeneId, geneId -> new LinkedList<>()).add(queryCigar);
                    }
                    var score = cigarSimilarity(targetCigar, queryCigar);

                    if (score > bestScore) {
                        bestScore = score;
                    }
                }
                firstRun = false;
            }

            if (bestScore >= SIMILARITY_CUTOFF) {
                results.add(genePair);
            } else {
                noSameTranscriptsLeftovers.add(genePair);

                var queryLeftOver = genePair.getQuery().getTranscripts().stream().map(transcript -> structureCigar(transcript, finalAllowedTypes)).collect(Collectors.toList());
                var targetLeftOver = genePair.getTarget().getTranscripts().stream().map(transcript -> structureCigar(transcript, finalAllowedTypes)).collect(Collectors.toList());

                queryLeftOverCigars.put(genePair.getQuery().getGeneId(), queryLeftOver);
                targetLeftOverCigars.put(genePair.getTarget().getGeneId(), targetLeftOver);

                unmappedTargetGenes.add(genePair.getTarget().getGeneId());
                unmappedQueryGenes.add(genePair.getQuery().getGeneId());
            }

        }

        var mappedNoSameIdWithCigar = new LinkedList<GenePair>();
        for (var leftOverQueryCigars : queryLeftOverCigars.entrySet()) {
            var queryGeneId = leftOverQueryCigars.getKey();
            var queryTranscriptCigars = leftOverQueryCigars.getValue();
            for (var leftOverTargetCigars : targetLeftOverCigars.entrySet()) {
                var targetGeneId = leftOverTargetCigars.getKey();
                var targetTranscriptCigars = leftOverTargetCigars.getValue();
                var bestScore = Double.NEGATIVE_INFINITY;

                if(targetGeneId.equals("ENSG00000162825") && queryGeneId.equals("ENSG00000122497")){
                    var a = 2;
                }

                for (var targetCigar : targetTranscriptCigars) {
                    for (var queryCigar : queryTranscriptCigars) {
                        var score = cigarSimilarity(targetCigar, queryCigar);
                        bestScore = Math.max(bestScore, score);
                    }
                }
                if (bestScore >= SIMILARITY_CUTOFF) {
                    mappedNoSameIdWithCigar.add(new GenePair(targetGtf.getGeneFeature(targetGeneId), queryGtf.getGeneFeature(queryGeneId)));

                    targetCigars.put(targetGeneId, targetTranscriptCigars);
                    queryCigars.put(queryGeneId, queryTranscriptCigars);

                    unmappedQueryGenes.remove(queryGeneId);
                    unmappedTargetGenes.remove(targetGeneId);
                }

            }
        }


        var queryPairsFoundComparingAgainstAllTargets = new LinkedList<GenePair>();
        var unmappedQueryCopy = List.copyOf(unmappedQueryGenes);
        for (var unmappedQueryGene : unmappedQueryCopy) {
            var queryGene = queryGtf.getGeneFeature(unmappedQueryGene);
            var unmappedQueryCigars = queryCigars.get(queryGene.getGeneId());

            if (unmappedQueryCigars == null) {
                unmappedQueryCigars = queryGene.getTranscripts().stream().map(t -> structureCigar(t, finalAllowedTypes)).toList();
                queryCigars.put(queryGene.getGeneId(), unmappedQueryCigars);
            }

            for (var unmappedQueryCigar : unmappedQueryCigars) {
                for (var targetCigarEntry : targetCigars.entrySet()) {
                    var targetTranscriptCigars = targetCigarEntry.getValue();
                    for (var targetCigar : targetTranscriptCigars) {

                        var score = cigarSimilarity(unmappedQueryCigar, targetCigar);

                        if (score >= SIMILARITY_CUTOFF) {
                            var targetGene = targetGtf.getGeneFeature(targetCigarEntry.getKey());
                            queryPairsFoundComparingAgainstAllTargets.add(new GenePair(targetGene, queryGene));

                            unmappedQueryGenes.remove(queryGene.getGeneId());
                            unmappedTargetGenes.remove(targetGene.getGeneId());

                            break;
                        }
                    }
                }
            }
        }


        var targetPairsFoundComparingAgainstAllQueries = new LinkedList<GenePair>();
        var unmappedTargetsCopy = List.copyOf(unmappedTargetGenes);
        for (var unmappedTargetGene : unmappedTargetsCopy) {
            var targetGene = targetGtf.getGeneFeature(unmappedTargetGene);
            var unmappedTargetCigars = targetCigars.get(targetGene.getGeneId());

            if (unmappedTargetCigars == null) {
                unmappedTargetCigars = targetGene.getTranscripts().stream().map(t -> structureCigar(t, finalAllowedTypes)).toList();
                targetCigars.put(targetGene.getGeneId(), unmappedTargetCigars);
            }

            for (var unmappedTargetCigar : unmappedTargetCigars) {
                for (var queryCigarEntry : queryCigars.entrySet()) {
                    var queryTranscriptCigars = queryCigarEntry.getValue();
                    for (var queryCigar : queryTranscriptCigars) {
                        if(targetGene.getGeneId().equals("ENSG00000162825") && queryCigarEntry.getKey().equals("ENSG00000122497")){
                            var a = 2;
                        }
                        var score = cigarSimilarity(queryCigar, unmappedTargetCigar);

                        if (score >= SIMILARITY_CUTOFF) {
                            var queryGene = queryGtf.getGeneFeature(queryCigarEntry.getKey());
                            targetPairsFoundComparingAgainstAllQueries.add(new GenePair(targetGene, queryGene));

                            unmappedQueryGenes.remove(queryGene.getGeneId());
                            unmappedTargetGenes.remove(targetGene.getGeneId());

                            break;
                        }
                    }
                }
            }
        }

        var numUnmapped = idMapping.size() - results.size() - mappedNoSameIdWithCigar.size();
        return List.of(results, mappedNoSameIdWithCigar, targetPairsFoundComparingAgainstAllQueries, queryPairsFoundComparingAgainstAllTargets);
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


}