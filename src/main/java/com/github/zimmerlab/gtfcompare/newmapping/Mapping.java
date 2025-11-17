package com.github.zimmerlab.gtfcompare.newmapping;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.zimmerlab.gtfcompare.model.GenePair;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class Mapping {
    public static List<GenePair> map(GtfFile targetGtf, GtfFile queryGtf) {
        var targetGeneIds = new HashSet<>(targetGtf.getAllGeneFeatureIds());
        var queryGeneIds = new HashSet<>(queryGtf.getAllGeneFeatureIds());

        var idMapping = idMapping(targetGtf, queryGtf, targetGeneIds, queryGeneIds);
        var noSameTranscripts = getGenesNoSameTranscript(idMapping);
        var leftOvers = getLeftoverGeneIds(targetGeneIds, queryGeneIds);


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

            for (var targetTranscript : target.getTranscripts()) {
                var targetTranscriptId = targetTranscript.getTranscriptId();
                var hasSameTranscript = false;
                for (var queryTranscript : query.getTranscripts()) {
                    var queryTranscriptId = queryTranscript.getTranscriptId();

                    if (queryTranscriptId.equals(targetTranscriptId)) {
                        hasSameTranscript = true;
                        break;
                    }
                }
                if (!hasSameTranscript) noSameTranscripts.add(genePair);
            }
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

}