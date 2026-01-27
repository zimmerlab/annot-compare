package com.github.zimmerlab.gtfcompare.newmappingval;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.zimmerlab.gtfcompare.model.GenePair;
import com.github.zimmerlab.gtfcompare.model.TranscriptPair;
import com.github.zimmerlab.gtfcompare.runner.NewMappingValidationRunner;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class MappingFileParser {
    private final List<TranscriptPair> mappings = new ArrayList<>();
    private final static Logger logger = LogManager.getLogger(MappingFileParser.class);
    public List<TranscriptPair> parse(String path, String contig, GtfFile queryGtfFile, GtfFile targetGtfFile) {
        try (BufferedReader reader = new BufferedReader(new FileReader(path))) {
            reader.readLine(); // HEADER
            reader.lines().forEach(line -> {
                parseLine(line, contig, queryGtfFile, targetGtfFile);
            });
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        return mappings;
    }

    private void parseLine(String line, String contig, GtfFile queryGtfFile, GtfFile targetGtfFile) {
        var parts = line.split("\t");
        var currentContig = parts[0];
        if(!currentContig.equals(contig)) return;
        var queryGeneId = parts[1];
        var targetGeneId = parts[2];
        var queryTranscriptId = parts[3];
        var targetTranscriptId = parts[4];

        var queryGene = queryGtfFile.getGeneFeature(queryGeneId);
        if (queryGene == null) {
            logger.warn("No gene found for query gene id {}", queryGeneId);
            return;
        };

        var targetGene = targetGtfFile.getGeneFeature(targetGeneId);
        if (targetGene == null) {
            logger.warn("No gene found for target gene id {}", targetGeneId);
            return;
        }

        var queryTranscriptObj = queryGene
                .getTranscripts().stream()
                .filter(transcript -> transcript.getTranscriptId().equals(queryTranscriptId)).toList();

        if(queryTranscriptObj.isEmpty()) {
            logger.warn("No transcript found for query transcript id {}. Gene Id {}", queryTranscriptId, queryGeneId);
            return;
        }

        var targetTranscriptObj = targetGene
                .getTranscripts().stream()
                .filter(transcript -> transcript.getTranscriptId().equals(targetTranscriptId)).toList();

        if(targetTranscriptObj.isEmpty()) {
            logger.warn("No transcript found for target transcript id {}. Gene Id {}", targetTranscriptId, targetGeneId);
            return;
        }

        mappings.add(new TranscriptPair(targetTranscriptObj.getFirst(), queryTranscriptObj.getFirst()));
    }

}
