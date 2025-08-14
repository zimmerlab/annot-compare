package com.github.zimmerlab.gtfcompare.utils;

import com.github.zimmerlab.gtfcompare.model.comparison.ComparisonResult;
import com.github.zimmerlab.gtfcompare.model.comparison.FeatureComparisonResult;
import com.github.zimmerlab.gtfcompare.model.comparison.RegionComparisonResult;
import com.github.zimmerlab.gtfcompare.model.comparison.TranscriptComparisonResult;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class ResultWriter {


    public static void writeComparisonResult(List<ComparisonResult> comparisonResults, String outputPath) {
        var lines = new ArrayList<OutputLine>();

        for (var result : comparisonResults) {
            collectGeneLines(lines, result);
        }

        var cmp = Comparator.comparingInt(OutputLine::getOrderKey).thenComparing(line -> String.join("\t", line.getColumns()));

        try (BufferedWriter writer = Files.newBufferedWriter(Path.of(outputPath), StandardOpenOption.CREATE, StandardOpenOption.APPEND)) {
            lines.stream().sorted(cmp).forEach(line -> {
                try {
                    writeLine(writer, line.getColumns());
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            });
        } catch (IOException e) {
            // TODO logging
            throw new RuntimeException(e);
        }
    }

    private static void collectGeneLines(List<OutputLine> lines, ComparisonResult result) {
        var targetGeneId = result.getTargetGeneId();
        var queryTargetId = result.getQueryGeneId();
        var geneComparisonResult = result.getGeneComparison();
        var queryBiotype = geneComparisonResult.getQueryBiotype();
        var targetBiotype = geneComparisonResult.getTargetBiotype();

        if(result.getTranscriptComparisons().size() > 1){
            var a = 2;
        }
        if (geneComparisonResult.isMissingInQueryFile()) {
            lines.add(new OutputLine(10, targetGeneId, queryTargetId,"","", "gene", "missingInQueryFile"));
            return;
        }

        if (geneComparisonResult.isMissingInTargetFile()) {
            lines.add(new OutputLine(10, targetGeneId, queryTargetId,"","", "gene", "missingInTargetFile"));
            return;
        }

        if (result.areSameGene()) {
            var tc = result.getTranscriptComparisons().get(0);
            String targetTranscriptId = "", queryTranscriptId = "";
            if(tc != null){
                targetTranscriptId = tc.getTargetTranscriptId();
                queryTranscriptId = tc.getQueryTranscriptId();
            }
            lines.add(new OutputLine(20, targetGeneId, queryTargetId, queryBiotype, targetBiotype, "", "", targetTranscriptId, queryTranscriptId));
        } else {

            if (geneComparisonResult.isStartDifferent())
                lines.add(new OutputLine(30, targetGeneId, queryTargetId, queryBiotype, targetBiotype, "gene", "start"));
            if (geneComparisonResult.isStopDifferent())
                lines.add(new OutputLine(40, targetGeneId, queryTargetId, queryBiotype, targetBiotype, "gene", "stop"));
            if (geneComparisonResult.isStrandDifferent())
                lines.add(new OutputLine(50, targetGeneId, queryTargetId, queryBiotype, targetBiotype, "gene", "strand"));
            if (!geneComparisonResult.getSequenceComparison().isSameSequence())
                lines.add(new OutputLine(60, targetGeneId, queryTargetId, queryBiotype, targetBiotype, "gene", "seq"));
            if (geneComparisonResult.isDifferentLength())
                lines.add(new OutputLine(70, targetGeneId, queryTargetId, queryBiotype, targetBiotype, "gene", "length"));
            if (geneComparisonResult.isContigDifferent())
                lines.add(new OutputLine(80, targetGeneId, queryTargetId, queryBiotype, targetBiotype, "gene", "contig"));

            for (var transcriptComparisonResult : result.getTranscriptComparisons()) {
                collectTranscriptLines(lines, targetGeneId, queryTargetId, transcriptComparisonResult, targetBiotype, queryBiotype);
            }
        }
    }

    private static void collectTranscriptLines(List<OutputLine> lines, String targetGeneId, String queryGeneId, TranscriptComparisonResult transcriptComparisonResult, String targetBiotype, String queryBiotype) {
        var targetTranscriptId = transcriptComparisonResult.getTargetTranscriptId();
        var queryTranscriptId = transcriptComparisonResult.getQueryTranscriptId();

        var queryTranscriptBiotype = transcriptComparisonResult.getQueryBiotype();
        var targetTranscriptBiotype = transcriptComparisonResult.getTargetBiotype();


        if (transcriptComparisonResult.isStartDifferent())
            lines.add(new OutputLine(100, targetGeneId, queryGeneId, targetBiotype, queryBiotype, "transcript", "start", targetTranscriptId, queryTranscriptId, targetTranscriptBiotype, queryTranscriptBiotype));
        if (transcriptComparisonResult.isStopDifferent())
            lines.add(new OutputLine(110, targetGeneId, queryGeneId, targetBiotype, queryBiotype, "transcript", "stop", targetTranscriptId, queryTranscriptId, targetTranscriptBiotype, queryTranscriptBiotype));
        if (transcriptComparisonResult.isSequenceDifferent())
            lines.add(new OutputLine(120, targetGeneId, queryGeneId, targetBiotype, queryBiotype, "transcript", "seq", targetTranscriptId, queryTranscriptId, targetTranscriptBiotype, queryTranscriptBiotype));
        if (transcriptComparisonResult.isTranscriptMissingInTargetGene())
            lines.add(new OutputLine(130, targetGeneId, queryGeneId, targetBiotype, queryBiotype, "transcript", "missingInFile1", targetTranscriptId, queryTranscriptId, targetTranscriptBiotype, queryTranscriptBiotype));
        if (transcriptComparisonResult.isTranscriptMissingInQueryGene())
            lines.add(new OutputLine(140, targetGeneId, queryGeneId, targetBiotype, queryBiotype, "transcript", "missingInFile2", targetTranscriptId, queryTranscriptId, targetTranscriptBiotype, queryTranscriptBiotype));
        if (transcriptComparisonResult.isLengthDifferent())
            lines.add(new OutputLine(150, targetGeneId, queryGeneId, targetBiotype, queryBiotype, "transcript", "length", targetTranscriptId, queryTranscriptId, targetTranscriptBiotype, queryTranscriptBiotype));

        for (FeatureComparisonResult fc : transcriptComparisonResult.getFeatureComparisons()) {
            collectFeatureLines(lines, targetGeneId, queryGeneId, targetTranscriptId, queryTranscriptId, fc, targetBiotype, queryBiotype, targetTranscriptBiotype, queryTranscriptBiotype);
        }
    }

    private static void collectFeatureLines(List<OutputLine> lines, String targetGeneId, String queryGeneId, String targetTranscriptId, String queryTranscriptId, FeatureComparisonResult featureComparisonResult, String targetBiotype, String queryBiotype, String targetTranscriptBiotype, String queryTranscriptBiotype) {
        var ft = featureComparisonResult.getFeatureType();

        if (featureComparisonResult.isMissingInTargetTranscript())
            lines.add(new OutputLine(200, targetGeneId, queryGeneId, targetBiotype, queryBiotype, ft, "missingInTranscript1", targetTranscriptId, queryTranscriptId, targetTranscriptBiotype, queryTranscriptBiotype));
        if (featureComparisonResult.isMissingInQueryTranscript())
            lines.add(new OutputLine(210, targetGeneId, queryGeneId, targetBiotype, queryBiotype, ft, "missingInTranscript2", targetTranscriptId, queryTranscriptId, targetTranscriptBiotype, queryTranscriptBiotype));

        for (RegionComparisonResult rc : featureComparisonResult.getRegionComparisons()) {
            if (rc.isLengthDifferenceFound())
                lines.add(new OutputLine(300, targetGeneId, queryGeneId, targetBiotype, queryBiotype, ft, "length", targetTranscriptId, queryTranscriptId, targetTranscriptBiotype, queryTranscriptBiotype));
            if (rc.isStartDifferent())
                lines.add(new OutputLine(310, targetGeneId, queryGeneId, targetBiotype, queryBiotype, ft, "start", targetTranscriptId, queryTranscriptId, targetTranscriptBiotype, queryTranscriptBiotype));
            if (rc.isEndDifferent())
                lines.add(new OutputLine(320, targetGeneId, queryGeneId, targetBiotype, queryBiotype, ft, "stop", targetTranscriptId, queryTranscriptId, targetTranscriptBiotype, queryTranscriptBiotype));
            if (rc.isSequenceDifferenceFound())
                lines.add(new OutputLine(330, targetGeneId, queryGeneId, targetBiotype, queryBiotype, ft, "seq", targetTranscriptId, queryTranscriptId, targetTranscriptBiotype, queryTranscriptBiotype));
            if (rc.isMissingInTargetFile())
                lines.add(new OutputLine(340, targetGeneId, queryGeneId, targetBiotype, queryBiotype, ft, "missingFeatureEntryFile1", targetTranscriptId, queryTranscriptId, targetTranscriptBiotype, queryTranscriptBiotype));
            if (rc.isMissingInQueryFile())
                lines.add(new OutputLine(350, targetGeneId, queryGeneId, targetBiotype, queryBiotype, ft, "missingFeatureEntryFile2", targetTranscriptId, queryTranscriptId, targetTranscriptBiotype, queryTranscriptBiotype));
        }
    }

    private static void writeLine(BufferedWriter writer, String... columns) throws IOException {
        writer.write(String.join("\t", columns));
        writer.newLine();
    }
}

