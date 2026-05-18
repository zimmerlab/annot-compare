package com.github.zimmerlab.gtfcompare.newmappingval;

import com.github.zimmerlab.gtfcompare.model.TranscriptPair;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

public class MappingValWriter {
    public static void write(String path, List<TranscriptPair> sameTranscripts, List<TranscriptPair> differentTranscripts){
        try(BufferedWriter writer = new BufferedWriter(new FileWriter(path))) {
            for (TranscriptPair pair : sameTranscripts) {
                writer.write(writeLine(pair, true));
            }

            for (TranscriptPair pair : differentTranscripts) {
                writer.write(writeLine(pair, false));
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private static String writeLine(TranscriptPair pair, boolean isSame){
        var targetTranscript = pair.getTarget();
        var queryTranscript = pair.getQuery();

        var targetGeneId = targetTranscript.getBaseData().getAttributes("gene_id").getFirst();
        var queryGeneId = queryTranscript.getBaseData().getAttributes("gene_id").getFirst();

        var targetTranscriptId = targetTranscript.getTranscriptId();
        var queryTranscriptId = queryTranscript.getTranscriptId();

        return String.format("%s\t%s\t%s\t%s\t%s\n",  targetGeneId, queryGeneId, targetTranscriptId, queryTranscriptId, isSame ? "true" : "false");
    }
}
