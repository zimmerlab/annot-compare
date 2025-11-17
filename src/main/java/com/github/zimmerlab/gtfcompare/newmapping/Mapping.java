package com.github.zimmerlab.gtfcompare.newmapping;

import com.github.kleinsamuel.gtfutils.GtfConfig;
import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.GeneFeature;
import com.github.zimmerlab.gtfcompare.model.GenePair;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class Mapping {
    public static List<GenePair> map(GtfFile targetGtf, GtfFile queryGtf) {
        var targetGeneIds = new HashSet<>(targetGtf.getAllGeneFeatureIds());
        var queryGeneIds = new HashSet<>(queryGtf.getAllGeneFeatureIds());

        var idMapping = idMapping(targetGtf, queryGtf, targetGeneIds, queryGeneIds);
        var noSameTranscripts = getGenesNoSameTranscript(idMapping);
        var leftOvers = getLeftoverGeneIds(targetGeneIds, queryGeneIds);

        for(var leftOver : leftOvers.targetOnly){
            var cigar = structureCigarString(targetGtf.getGeneFeature(leftOver));
        }

        for(var leftOver : leftOvers.queryOnly){
            var cigar = structureCigarString(queryGtf.getGeneFeature(leftOver));
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

    private record LeftoverGeneIds(Set<String> targetOnly, Set<String> queryOnly){}
    private static LeftoverGeneIds getLeftoverGeneIds(Set<String> targetGeneIds, Set<String> queryGeneIds){
        var targetOnly = new HashSet<>(targetGeneIds);
        targetOnly.removeAll(queryGeneIds);

        var queryOnly = new HashSet<>(queryGeneIds);
        queryOnly.removeAll(targetGeneIds);

        return new LeftoverGeneIds(targetOnly, queryOnly);
    }

    private static String structureCigarString(GeneFeature geneFeature){
        var geneExons = geneFeature.getTranscripts()
                .stream()
                .flatMap(t -> t.getFeatures(GtfConfig.TYPE_EXON_DEFAULT).stream())
                .collect(Collectors.toSet());

        StringBuilder sb = new StringBuilder();
        var lastStop = 0;
        var isFirst = true;
        for(var geneExon : geneExons){
            var baseData = geneExon.getBaseData();
            var start = baseData.getStart();
            var stop = baseData.getEnd();

            if(!isFirst){
                sb.append("I");
                sb.append((start - lastStop + 1));
            }

            var length = stop - start + 1;

            sb.append("E");
            sb.append(length);

            isFirst = false;
            lastStop = stop;
        }

        return sb.toString();
    }
}