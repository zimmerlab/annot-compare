package com.github.zimmerlab.gtfcompare.analysis.lifted;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.zimmerlab.gtfcompare.model.GenePair;
import com.github.zimmerlab.gtfcompare.utils.GenomeSequenceExtractor;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class LiftedDifferences {
    private final static Logger logger = LogManager.getLogger(LiftedDifferences.class);
    private static boolean headerAdded = false;
    private static class Stats {
        int numQueryGenes;
        int numTargetGenes;
        int numLiftedQueryGenes;
        int normalMappingSize;
        int liftedMappingSize;
        int numDifferentPairsByPosition;
        int numDifferentPairsBySequence;
        List<String> notFoundInLifted = new ArrayList<>();
    }

    private static final Map<String, Stats> results = new HashMap<>();

    public static void analyze(String contig, GtfFile queryGtf, GtfFile targetGtf, GtfFile liftedQueryGtf, GenomeSequenceExtractor queryGSE, GenomeSequenceExtractor targetGSE, String outputPath) {
        try (BufferedWriter byPositionWriter = new BufferedWriter(new FileWriter(outputPath + ".diffByPosition"));
             BufferedWriter bySeqWriter = new BufferedWriter(new FileWriter(outputPath + ".diffBySeq"));
             BufferedWriter normalMappingWriter = new BufferedWriter(new FileWriter(outputPath + ".normalMapping"));
             BufferedWriter liftedMappingWriter = new BufferedWriter(new FileWriter(outputPath + ".liftedMapping"));
             BufferedWriter missingInLiftedWriter = new BufferedWriter(new FileWriter(outputPath + ".missingInLifted"));) {

        addHeader(byPositionWriter, bySeqWriter, normalMappingWriter, liftedMappingWriter, missingInLiftedWriter);

        var s = new Stats();

        s.numQueryGenes = queryGtf.getAllGeneFeatureIds().size();
        s.numTargetGenes = targetGtf.getAllGeneFeatureIds().size();
        s.numLiftedQueryGenes = liftedQueryGtf.getAllGeneFeatureIds().size();

        var mappingWithNormal = find1To1Mapping(queryGtf, targetGtf);
        writeGenePairs(normalMappingWriter, mappingWithNormal);

        var mappingWithLifted = find1To1Mapping(liftedQueryGtf, targetGtf);
        writeGenePairs(liftedMappingWriter, mappingWithLifted);

        s.normalMappingSize = mappingWithNormal.size();
        s.liftedMappingSize = mappingWithLifted.size();

        var differentPairsByPosition = findDifferentPairsByPosition(mappingWithLifted);
        writeGenePairs(byPositionWriter, differentPairsByPosition);
        s.numDifferentPairsByPosition = differentPairsByPosition.size();

        var differentPairsBySequence = findDifferentPairsBySequence(mappingWithLifted, queryGtf, queryGSE, targetGSE, s);
        writeGenePairs(bySeqWriter, differentPairsBySequence);
        s.numDifferentPairsBySequence = differentPairsBySequence.size();

        results.put(contig, s);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private static void addHeader(BufferedWriter... writers) throws IOException {
        if(!headerAdded){
            var headerString = "contig\tqueryGeneId\ttargetGeneId\n";
            for(var writer : writers){
                writer.write(headerString);
            }
            headerAdded = true;
        }
    }

    public static void printResults(String outputPath) {

        try (BufferedWriter overviewWriter = new BufferedWriter(new FileWriter(outputPath));
             BufferedWriter missingInLiftedWriter = new BufferedWriter(new FileWriter(outputPath + ".missingInLifted"))) {

            overviewWriter.write("Chrom\tQueryGenes\tTargetGenes\tLiftedQueryGenes\t1to1Normal\t1to1Lifted\tDiffPos\tDiffSeq\n");

            int totalQ = 0, totalT = 0, totalN = 0, totalL = 0, totalP = 0, totalS = 0;

            for (var entry : results.entrySet()) {
                String chr = entry.getKey();
                Stats s = entry.getValue();

                overviewWriter.write(chr + "\t" + s.numQueryGenes + "\t" + s.numTargetGenes + "\t" + s.numLiftedQueryGenes + "\t" + s.normalMappingSize + "\t" + s.liftedMappingSize + "\t" + s.numDifferentPairsByPosition + "\t" + s.numDifferentPairsBySequence + "\n");

                totalQ += s.numQueryGenes;
                totalT += s.numTargetGenes;
                totalN += s.normalMappingSize;
                totalL += s.liftedMappingSize;
                totalP += s.numDifferentPairsByPosition;
                totalS += s.numDifferentPairsBySequence;

                writeGenes(missingInLiftedWriter, s.notFoundInLifted);
            }

            overviewWriter.write("TOTAL\t" + totalQ + "\t" + totalT + "\t" + totalN + "\t" + totalL + "\t" + totalP + "\t" + totalS + "\n");
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private static void writeGenePairs(BufferedWriter fw, List<GenePair> genePairs) throws IOException {
        for (var genePair : genePairs) {
            var targetGeneId = genePair.getTarget().getGeneId();
            var queryGeneId = genePair.getQuery().getGeneId();
            var contig = genePair.getTarget().getBaseData().getContig();

            fw.write(contig + "\t" + queryGeneId + "\t" + targetGeneId + "\n");
        }
    }

    private static void writeGenes(BufferedWriter fw, List<String> genes) throws IOException {
        for (var gene : genes) {
            fw.write(gene + "\n");
        }
    }

    private static List<GenePair> find1To1Mapping(GtfFile queryGtf, GtfFile targetGtf) {
        var mapping = new ArrayList<GenePair>();

        for (var geneId : queryGtf.getAllGeneFeatureIds()) {
            var targetGene = targetGtf.getGeneFeature(geneId);
            if (targetGene == null) continue;

            mapping.add(new GenePair(targetGene, queryGtf.getGeneFeature(geneId)));
        }

        return mapping;
    }

    private static List<GenePair> findDifferentPairsByPosition(List<GenePair> mapping) {
        var differentPairs = new ArrayList<GenePair>();
        for (var genePair : mapping) {
            var targetGeneBaseData = genePair.getTarget().getBaseData();
            var queryGeneBaseData = genePair.getQuery().getBaseData();

            if (targetGeneBaseData.getStart() != queryGeneBaseData.getStart()) {
                differentPairs.add(genePair);
                continue;
            }

            if (targetGeneBaseData.getEnd() != queryGeneBaseData.getEnd()) {
                differentPairs.add(genePair);
            }
        }

        return differentPairs;
    }

    private static List<GenePair> findDifferentPairsBySequence(List<GenePair> mapping, GtfFile queryGtf, GenomeSequenceExtractor queryGSE, GenomeSequenceExtractor targetGSE, Stats s) {
        var differentPairs = new ArrayList<GenePair>();
        for (var genePair : mapping) {
            var targetGeneBaseData = genePair.getTarget().getBaseData();
            var queryGeneBaseData = genePair.getQuery().getBaseData();

            if (targetGeneBaseData.getStart() != queryGeneBaseData.getStart()) {
                continue;
            }

            if (targetGeneBaseData.getEnd() != queryGeneBaseData.getEnd()) {
                continue;
            }

            var queryGtfGeneId = genePair.getQuery().getGeneId();
            var queryGtfGene = queryGtf.getGeneFeature(queryGtfGeneId);
            if (queryGtfGene == null) {
                logger.debug("{} not found", queryGtfGeneId);
                s.notFoundInLifted.add(queryGtfGeneId);
                continue;
            }

            var queryStart = queryGtfGene.getBaseData().getStart();
            var queryEnd = queryGtfGene.getBaseData().getEnd();

            try {
                var targetSeq = targetGSE.getSequence(targetGeneBaseData.getContig(), targetGeneBaseData.getStart(), targetGeneBaseData.getEnd());
                var querySeq = queryGSE.getSequence(queryGtfGene.getBaseData().getContig(), queryStart, queryEnd);

                if (targetSeq.equals(querySeq)) continue;

                differentPairs.add(genePair);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }

        return differentPairs;
    }
}
