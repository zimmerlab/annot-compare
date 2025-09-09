package com.github.zimmerlab.gtfcompare.mapping;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.GeneFeature;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.model.GenePair;
import com.github.zimmerlab.gtfcompare.model.config.ComparisonConfig;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class MapById {
    private List<GenePair> getGenePairsByExactId(GtfFile targetGtf, GtfFile queryGtf, ComparisonConfig config) {
        var targetGeneMap = new HashMap<String, GeneFeature>();
        var queryGeneMap = new HashMap<String, GeneFeature>();

        var everyBiotypeAllowed = config.getAllowedGeneBiotypes().isEmpty();
        for (var geneId : targetGtf.getAllGeneFeatureIds()) {
            var gene = targetGtf.getGeneFeature(geneId);

            if (everyBiotypeAllowed) {
                targetGeneMap.put(geneId, gene);
                continue;
            }
            var baseData = gene.getBaseData();
            var geneBiotypeEmpty = baseData.getAttributes("gene_biotype").isEmpty();


            String geneBiotype = geneBiotypeEmpty ? baseData.getSource() : baseData.getAttributes("gene_biotype").getFirst();
            if (config.getAllowedGeneBiotypes().contains(geneBiotype)) targetGeneMap.put(geneId, gene);
        }

        for (var geneId : queryGtf.getAllGeneFeatureIds()) {
            var gene = queryGtf.getGeneFeature(geneId);
            if (everyBiotypeAllowed) {
                queryGeneMap.put(geneId, gene);
                continue;
            }
            var baseData = gene.getBaseData();
            var geneBiotypeEmpty = baseData.getAttributes("gene_biotype").isEmpty();

            String geneBiotype = geneBiotypeEmpty ? baseData.getSource() : baseData.getAttributes("gene_biotype").getFirst();
            if (config.getAllowedGeneBiotypes().contains(geneBiotype)) queryGeneMap.put(geneId, gene);
        }

        var allGeneIds = new HashSet<String>();
        allGeneIds.addAll(targetGeneMap.keySet());
        allGeneIds.addAll(queryGeneMap.keySet());

        var genePairs = new ArrayList<GenePair>();
        for (var geneId : allGeneIds) {
            GeneFeature t1 = targetGeneMap.get(geneId);
            GeneFeature t2 = queryGeneMap.get(geneId);
            genePairs.add(new GenePair(t1, t2));
        }
        return genePairs;
    }

    public static List<TranscriptFeature> changeTargetAttributesBasedOnLiftoff(GtfFile targetGtf, GtfFile liftOffGtf) {
        var targetTranscriptMap = new HashMap<String, TranscriptFeature>();
        var liftOffTranscriptMap = new HashMap<String, TranscriptFeature>();


        for (var geneId : liftOffGtf.getAllGeneFeatureIds()) {
            var gene = liftOffGtf.getGeneFeature(geneId);
            for(var transcript : gene.getTranscripts()) {
                liftOffTranscriptMap.put(transcript.getTranscriptId(), transcript);
            }
        }

        for (var geneId : targetGtf.getAllGeneFeatureIds()) {
            var gene = targetGtf.getGeneFeature(geneId);
            for(var transcript : gene.getTranscripts()) {
                targetTranscriptMap.put(transcript.getTranscriptId(), transcript);
            }
        }

        var allTranscriptIds = new HashSet<String>();
        allTranscriptIds.addAll(targetTranscriptMap.keySet());
        allTranscriptIds.addAll(liftOffTranscriptMap.keySet());

        var transcriptPairs = new ArrayList<TranscriptFeature>();
        for (var transcriptIds : allTranscriptIds) {
            var t1 = targetTranscriptMap.get(transcriptIds);
            var t2 = liftOffTranscriptMap.get(transcriptIds);

            var targetBaseData = t1.getBaseData();
            var liftOffBaseData = t2.getBaseData();

            targetBaseData.setForwardStrand(liftOffBaseData.isForwardStrand());
            // change necessary attributes... method prob. not necessary, but you get the gist

            transcriptPairs.add(t1);
        }
        return transcriptPairs;
    }
}
