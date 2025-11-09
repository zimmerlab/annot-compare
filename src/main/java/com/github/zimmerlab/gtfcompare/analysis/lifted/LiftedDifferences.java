package com.github.zimmerlab.gtfcompare.analysis.lifted;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.zimmerlab.gtfcompare.model.GenePair;

import java.util.ArrayList;
import java.util.List;

public class LiftedDifferences {

    private static int numQueryGenes;
    private static int numTargetGenes;
    private static int normalMappingSize;
    private static int liftedMappingSize;
    private static int numDifferentPairsByPosition;

    public static void analyze(GtfFile queryGtf, GtfFile targetGtf, GtfFile liftedQueryGtf){
        numQueryGenes = queryGtf.getAllGeneFeatureIds().size();
        numTargetGenes = targetGtf.getAllGeneFeatureIds().size();

        var mappingWithNormal = find1To1Mapping(queryGtf, targetGtf);
        var mappingWithLifted = find1To1Mapping(liftedQueryGtf, targetGtf);

        normalMappingSize = mappingWithNormal.size();
        liftedMappingSize = mappingWithLifted.size();

        var differentPairsByPosition = findDifferentPairsByPosition(mappingWithLifted);
        numDifferentPairsByPosition = differentPairsByPosition.size();
    }

    private static List<GenePair> find1To1Mapping(GtfFile queryGtf, GtfFile targetGtf){
        var mapping = new ArrayList<GenePair>();

        for(var geneId : queryGtf.getAllGeneFeatureIds()){
            var targetGene = targetGtf.getGeneFeature(geneId);
            if(targetGene == null) continue;

            mapping.add(new GenePair(targetGene, queryGtf.getGeneFeature(geneId)));
        }

        return mapping;
    }

    private static List<GenePair> findDifferentPairsByPosition(List<GenePair> mapping){
        var differentPairs = new ArrayList<GenePair>();
        for(var genePair : mapping){
            var targetGeneBaseData = genePair.getTargetGene().getBaseData();
            var queryGeneBaseData = genePair.getQueryGene().getBaseData();

            if(targetGeneBaseData.getStart() != queryGeneBaseData.getStart()) {
                differentPairs.add(genePair);
                continue;
            }

            if(targetGeneBaseData.getEnd() != queryGeneBaseData.getEnd()) {
                differentPairs.add(genePair);
            }
        }

        return differentPairs;
    }
}
