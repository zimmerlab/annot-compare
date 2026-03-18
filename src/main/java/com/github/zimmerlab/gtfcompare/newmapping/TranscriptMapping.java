package com.github.zimmerlab.gtfcompare.newmapping;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.model.TranscriptPair;
import com.github.zimmerlab.gtfcompare.newmapping.model.*;
import com.github.zimmerlab.gtfcompare.utils.GenomeSequenceExtractor;

import java.util.*;
import java.util.function.Consumer;
import java.util.stream.Collectors;

public class TranscriptMapping {


    public static UnmappedResult map(List<ResultWithOrigin> mappings, GtfFile targetGtf, GtfFile queryGtf, Set<String> allowedTypes, GenomeSequenceExtractor targetExtractor, GenomeSequenceExtractor queryExtractor, Consumer<ResultWithOrigin> resultConsumer) {
        var unmappedQueriesWithTranscriptNames = new HashMap<String, String>();
        var unmappedTargetsWithTranscriptNames = new HashMap<String, String>();
        for (var mapping : mappings) {

            var genePair = mapping.genePair();

            var targetGene = genePair.getTarget();
            var queryGene = genePair.getQuery();

            var transcriptMappings = mapTranscripts(targetGene.getTranscripts(), queryGene.getTranscripts(), allowedTypes, targetExtractor, queryExtractor);



            var currentUnmappedQueriesWithTranscriptNames = Utils.addTranscriptNames(queryGene, transcriptMappings.unmappedQueries());
            var currentUnmappedTargetsWithTranscriptNames = Utils.addTranscriptNames(targetGene, transcriptMappings.unmappedTargets());

            unmappedQueriesWithTranscriptNames.putAll(currentUnmappedQueriesWithTranscriptNames);
            unmappedTargetsWithTranscriptNames.putAll(currentUnmappedTargetsWithTranscriptNames);

            for (var transcriptMapping : transcriptMappings.matches()) {
                var targetTranscriptId = transcriptMapping.target();
                var queryTranscriptId = transcriptMapping.query();
                var currentOrigins = new ArrayList<>(mapping.origins());
                if(targetTranscriptId.equals(queryTranscriptId)) {
                    currentOrigins.add(MappingOrigin.TRANSCRIPT_ID_MAPPING);
                }

                var resultWithOrigin = new ResultWithOrigin(genePair, currentOrigins, mapping.geneDistance(), targetTranscriptId, queryTranscriptId);
                resultConsumer.accept(resultWithOrigin);
            }
        }

        return new UnmappedResult(unmappedQueriesWithTranscriptNames, unmappedTargetsWithTranscriptNames);
    }

    private static TranscriptMappingResult mapTranscripts(List<TranscriptFeature> targetTranscripts, List<TranscriptFeature> queryTranscripts, Set<String> allowedTypes, GenomeSequenceExtractor targetExtractor, GenomeSequenceExtractor queryExtractor) {
        var matches = new ArrayList<TranscriptIdPair>();
        var queryTranscriptCache = new HashMap<String, List<CigarOp>>();

        var matchedTargetIds = new HashSet<String>();
        var matchedQueryIds = new HashSet<String>();

        var allTargetIds = targetTranscripts.stream().map(TranscriptFeature::getTranscriptId).collect(Collectors.toSet());
        var allQueryIds = queryTranscripts.stream().map(TranscriptFeature::getTranscriptId).collect(Collectors.toSet());

        for (var targetTranscript : targetTranscripts) {
            var targetCigar = Similarity.structureCigar(targetTranscript, allowedTypes, targetExtractor);
            var targetTranscriptId = targetTranscript.getTranscriptId();

            for (var queryTranscript : queryTranscripts) {
                var queryTranscriptId = queryTranscript.getTranscriptId();

                var queryCigar = queryTranscriptCache.computeIfAbsent(queryTranscriptId, k -> Similarity.structureCigar(queryTranscript, allowedTypes, queryExtractor));

                if (!Similarity.isSimilar(targetCigar, queryCigar)) continue;

                matches.add(new TranscriptIdPair(targetTranscriptId, queryTranscriptId));

                matchedTargetIds.add(targetTranscriptId);
                matchedQueryIds.add(queryTranscriptId);
            }
        }

        var unmappedTargets = allTargetIds.stream().filter(id -> !matchedTargetIds.contains(id)).toList();

        var unmappedQueries = allQueryIds.stream().filter(id -> !matchedQueryIds.contains(id)).toList();

        return new TranscriptMappingResult(matches, unmappedTargets, unmappedQueries);
    }
}
