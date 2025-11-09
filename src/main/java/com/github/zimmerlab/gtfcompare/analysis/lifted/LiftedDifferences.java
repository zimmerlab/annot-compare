package com.github.zimmerlab.gtfcompare.analysis.lifted;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.zimmerlab.gtfcompare.model.GenePair;
import com.github.zimmerlab.gtfcompare.utils.GenomeSequenceExtractor;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class LiftedDifferences {

    private static int numQueryGenes;
    private static int numTargetGenes;
    private static int normalMappingSize;
    private static int liftedMappingSize;
    private static int numDifferentPairsByPosition;
    private static int numDifferentPairsBySequence;

    public static void analyze(GtfFile queryGtf, GtfFile targetGtf, GtfFile liftedQueryGtf, GenomeSequenceExtractor queryGSE, GenomeSequenceExtractor targetGSE){
        numQueryGenes = queryGtf.getAllGeneFeatureIds().size();
        numTargetGenes = targetGtf.getAllGeneFeatureIds().size();

        var mappingWithNormal = find1To1Mapping(queryGtf, targetGtf);
        var mappingWithLifted = find1To1Mapping(liftedQueryGtf, targetGtf);

        normalMappingSize = mappingWithNormal.size();
        liftedMappingSize = mappingWithLifted.size();

        var differentPairsByPosition = findDifferentPairsByPosition(mappingWithLifted);
        numDifferentPairsByPosition = differentPairsByPosition.size();

        findDifferentPairsBySequence(mappingWithLifted, queryGtf, queryGSE, targetGSE);
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

    private static List<GenePair> findDifferentPairsBySequence(List<GenePair> mapping, GtfFile queryGtf, GenomeSequenceExtractor queryGSE, GenomeSequenceExtractor targetGSE){
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
                continue;
            }

            var queryGtfGene = queryGtf.getGeneFeature(genePair.getQueryGene().getGeneId());
            if(queryGtfGene == null) {
                differentPairs.add(genePair);
                System.out.println("wtf");
                continue;
            }

            var queryStart = queryGtfGene.getBaseData().getStart();
            var queryEnd = queryGtfGene.getBaseData().getEnd();

            try {
                var targetSeq = targetGSE.getSequence(targetGeneBaseData.getContig(), targetGeneBaseData.getStart(), targetGeneBaseData.getEnd());
                var querySeq = queryGSE.getSequence(queryGtfGene.getBaseData().getContig(), queryStart, queryEnd);

                if(targetSeq.equals(querySeq)) continue;

                differentPairs.add(genePair);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }

        return differentPairs;
    }
}
