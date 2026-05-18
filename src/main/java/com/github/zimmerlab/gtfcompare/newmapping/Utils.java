package com.github.zimmerlab.gtfcompare.newmapping;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.GeneFeature;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;

public class Utils {

    public static String getAttribute(GeneFeature gene, String key) {
        var attrs = gene.getBaseData().getAttributes(key);

        if (attrs == null || attrs.isEmpty()) {
            return null;
        }

        return attrs.getFirst();
    }


    public static String getName(GtfFile gtfFile, String geneId) {
        var gene = gtfFile.getGeneFeature(geneId);
        var name = getAttribute(gene, "gene_name");
        return name != null ? name : "";
    }

    public static Map<String, String> addNames(Collection<String> geneIds, GtfFile gtfFile) {
        var validGeneIds = new ArrayList<String>();
        for (String geneId : geneIds) {
            var gene = gtfFile.getGeneFeature(geneId);

            var geneBiotype = getAttribute(gene, "gene_biotype");
            var source = gene.getBaseData().getSource();

            if (geneBiotype == null || (!geneBiotype.equals("protein_coding") && !source.equals("protein_coding")))
                continue;

            validGeneIds.add(geneId);
        }

        return validGeneIds.stream()
                .collect(Collectors.toMap(
                        id -> id,
                        id -> getName(gtfFile, id)
                ));
    }

    public static Map<String, String> addTranscriptNames(GeneFeature gene, Collection<String> transcriptIds) {
        var resultMap = new HashMap<String, String>();
        for (var transcript : gene.getTranscripts()) {
            var transcriptId = transcript.getTranscriptId();
            if (!transcriptIds.contains(transcriptId)) continue;

            var transcriptNameList = transcript.getBaseData().getAttributes("transcript_name");
            if(transcriptNameList == null || transcriptNameList.isEmpty()){
                resultMap.put(transcriptId, "");
            } else{
                resultMap.put(transcriptId, transcriptNameList.getFirst());
            }
        }

        return resultMap;
    }

}
