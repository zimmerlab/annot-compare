package com.github.zimmerlab.gtfcompare.newmapping;

import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.newmapping.model.*;
import com.github.zimmerlab.gtfcompare.utils.GenomeSequenceExtractor;

import java.util.*;
import java.util.function.Consumer;
import java.util.stream.Collectors;

public class TranscriptMapping {


    public static UnmappedResult map(List<ResultWithOrigin> mappings, Set<String> allowedTypes, GenomeSequenceExtractor targetExtractor, GenomeSequenceExtractor queryExtractor, Consumer<ResultWithOrigin> resultConsumer) {
        var mappedQueryTranscriptIds = new HashSet<String>();
        var mappedTargetTranscriptIds = new HashSet<String>();

        var seenQueryTranscripts = new HashMap<String, String>();
        var seenTargetTranscripts = new HashMap<String, String>();

        for (var mapping : mappings) {
            var genePair = mapping.genePair();

            var targetGene = genePair.getTarget();
            var queryGene = genePair.getQuery();

            seenQueryTranscripts.putAll(Utils.addTranscriptNames(queryGene, queryGene.getTranscripts().stream().map(TranscriptFeature::getTranscriptId).collect(Collectors.toSet())));
            seenTargetTranscripts.putAll(Utils.addTranscriptNames(targetGene, targetGene.getTranscripts().stream().map(TranscriptFeature::getTranscriptId).collect(Collectors.toSet())));

            var transcriptMappings = mapTranscripts(targetGene.getTranscripts(), queryGene.getTranscripts(), allowedTypes, targetExtractor, queryExtractor);

            for (var transcriptMapping : transcriptMappings) {
                var targetTranscriptId = transcriptMapping.target();
                var queryTranscriptId = transcriptMapping.query();

                mappedTargetTranscriptIds.add(targetTranscriptId);
                mappedQueryTranscriptIds.add(queryTranscriptId);

                var currentOrigins = new ArrayList<>(mapping.origins());
                if (targetTranscriptId.equals(queryTranscriptId)) {
                    currentOrigins.add(MappingOrigin.TRANSCRIPT_ID_MAPPING);
                }

                var resultWithOrigin = new ResultWithOrigin(genePair, currentOrigins, mapping.geneDistance(), targetTranscriptId, queryTranscriptId);
                resultConsumer.accept(resultWithOrigin);
            }
        }

        var unmappedQueriesWithTranscriptNames = new HashMap<String, String>();
        for (var entry : seenQueryTranscripts.entrySet()) {
            if (!mappedQueryTranscriptIds.contains(entry.getKey())) {
                unmappedQueriesWithTranscriptNames.put(entry.getKey(), entry.getValue());
            }
        }

        var unmappedTargetsWithTranscriptNames = new HashMap<String, String>();
        for (var entry : seenTargetTranscripts.entrySet()) {
            if (!mappedTargetTranscriptIds.contains(entry.getKey())) {
                unmappedTargetsWithTranscriptNames.put(entry.getKey(), entry.getValue());
            }
        }

        return new UnmappedResult(unmappedQueriesWithTranscriptNames, unmappedTargetsWithTranscriptNames);
    }

    private static List<TranscriptIdPair> mapTranscripts(List<TranscriptFeature> targetTranscripts, List<TranscriptFeature> queryTranscripts, Set<String> allowedTypes, GenomeSequenceExtractor targetExtractor, GenomeSequenceExtractor queryExtractor) {
        var matches = new ArrayList<TranscriptIdPair>();
        var queryTranscriptCache = new HashMap<String, List<CigarOp>>();

        for (var targetTranscript : targetTranscripts) {
            var targetCigar = Similarity.structureCigar(targetTranscript, allowedTypes, targetExtractor);
            var targetTranscriptId = targetTranscript.getTranscriptId();

            for (var queryTranscript : queryTranscripts) {
                var queryTranscriptId = queryTranscript.getTranscriptId();

                var queryCigar = queryTranscriptCache.computeIfAbsent(queryTranscriptId, k -> Similarity.structureCigar(queryTranscript, allowedTypes, queryExtractor));

                if (!Similarity.isSimilar(targetCigar, queryCigar)) continue;

                matches.add(new TranscriptIdPair(targetTranscriptId, queryTranscriptId));
            }
        }

        return matches;
    }
}
