package com.github.zimmerlab.gtfcompare.newmapping;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
import com.github.zimmerlab.gtfcompare.model.GenePair;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

public class IdMapping {
    public static List<GenePair> map(GtfFile targetGtf, GtfFile queryGtf) {
        var targetGeneIds = new HashSet<>(targetGtf.getAllGeneFeatureIds());
        var queryGeneIds = new HashSet<>(queryGtf.getAllGeneFeatureIds());

        targetGeneIds.retainAll(queryGeneIds);
        var mapping = new ArrayList<GenePair>(targetGeneIds.size());

        for(var geneId : targetGeneIds){
            mapping.add(new GenePair(targetGtf.getGeneFeature(geneId), queryGtf.getGeneFeature(geneId)));
        }
        return mapping;
    }
}
