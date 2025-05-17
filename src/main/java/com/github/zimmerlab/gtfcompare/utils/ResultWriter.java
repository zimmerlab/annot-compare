package com.github.zimmerlab.gtfcompare.utils;

import com.github.zimmerlab.gtfcompare.model.comparison.ComparisonResult;
import com.github.zimmerlab.gtfcompare.model.comparison.FeatureComparisonResult;
import com.github.zimmerlab.gtfcompare.model.comparison.RegionComparisonResult;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

public class ResultWriter {

    public static void writeComparisonResult(List<ComparisonResult> comparisonResults, String outputPath) {
        try (var writer = new BufferedWriter(new FileWriter((outputPath)))) {
            writer.write("targetGeneId\tqueryGeneId\tcategory\tdifference\n");
            for (var comparisonResult : comparisonResults) {

                var targetGeneId = comparisonResult.getTargetGeneId();
                var queryGeneId = comparisonResult.getQueryGeneId();

                var geneComparison = comparisonResult.getGeneComparison();

                if (geneComparison.isMissingInQueryFile()) {
                    writer.write(targetGeneId + "\t" + queryGeneId + "\tgene\tmissingInQueryFile\n");
                    continue;
                }



                if (comparisonResult.areSameGene()) {
                    writer.write(targetGeneId + "\t" + queryGeneId + "\n");
                } else {

                    if (geneComparison.isStartDifferent()) {
                        writer.write(targetGeneId + "\t" + queryGeneId + "\tgene\tstart\n");
                    }
                    if (geneComparison.isStopDifferent()) {
                        writer.write(targetGeneId + "\t" + queryGeneId + "\tgene\tstop\n");
                    }
                    if (geneComparison.isStrandDifferent()) {
                        writer.write(targetGeneId + "\t" + queryGeneId + "\tgene\tstrand\n");
                    }
                    if (!geneComparison.getSequenceComparison().isSameSequence()) {
                        writer.write(targetGeneId + "\t" + queryGeneId+ "\tgene\tseq\n");
                    }
                    if (geneComparison.isDifferentLength()) {
                        writer.write(targetGeneId + "\t" + queryGeneId+ "\tgene\tlength\n");
                    }

                    var transcriptComparisonResult = comparisonResult.getTranscriptComparisons();

                    for (var transcriptComparison : transcriptComparisonResult) {
                        if (transcriptComparison.isStartDifferent()) {
                            writer.write(targetGeneId + "\t" + queryGeneId+ "\ttranscript\tstart_" + transcriptComparison.getQueryTranscriptId() + "\n");
                        }
                        if (transcriptComparison.isStopDifferent()) {
                            writer.write(targetGeneId + "\t" + queryGeneId+ "\ttranscript\tstop_" + transcriptComparison.getQueryTranscriptId() + "\n");
                        }
                        if (transcriptComparison.isSequenceDifferent()) {
                            writer.write(targetGeneId + "\t" + queryGeneId+ "\ttranscript\tseq_" + transcriptComparison.getQueryTranscriptId() + "\n");
                        }
                        if (transcriptComparison.isTranscriptMissingInTargetGene()) {
                            writer.write(targetGeneId + "\t" + queryGeneId+ "\ttranscript\tmissingInFile1_" + transcriptComparison.getQueryTranscriptId() + "\n");
                        }
                        if (transcriptComparison.isTranscriptMissingInQueryGene()) {
                            writer.write(targetGeneId + "\t" + queryGeneId+ "\ttranscript\tmissingInFile2_" + transcriptComparison.getQueryTranscriptId() + "\n");
                        }
                        if (transcriptComparison.isLengthDifferent()) {
                            writer.write(targetGeneId + "\t" + queryGeneId+ "\ttranscript\tlength_" + transcriptComparison.getQueryTranscriptId() + "\n");

                        }

                        for (FeatureComparisonResult featureComparison : transcriptComparison.getFeatureComparisons()) {
                            if (featureComparison.isMissingInTargetTranscript()) {
                                writer.write(targetGeneId + "\t" + queryGeneId+ "\tfeature\tmissingInTranscript1_" + transcriptComparison.getQueryTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                            }
                            if (featureComparison.isMissingInQueryTranscript()) {
                                writer.write(targetGeneId + "\t" + queryGeneId+ "\tfeature\tmissingInTranscript2_" + transcriptComparison.getQueryTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                            }

                            for (RegionComparisonResult regionComparison : featureComparison.getRegionComparisons()) {
                                if (regionComparison.isLengthDifferenceFound()) {
                                    writer.write(targetGeneId + "\t" + queryGeneId+ "\tfeature\tlength_" + transcriptComparison.getQueryTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                }
                                if (regionComparison.isPositionDifferenceFound()) {
                                    writer.write(targetGeneId + "\t" + queryGeneId+ "\tfeature\tposition_" + transcriptComparison.getQueryTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                }
                                if (regionComparison.isSequenceDifferenceFound()) {
                                    writer.write(targetGeneId + "\t" + queryGeneId+ "\tfeature\tseq_" + transcriptComparison.getQueryTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                }
                                if (regionComparison.isMissingInTargetFile()) {
                                    writer.write(targetGeneId + "\t" + queryGeneId+ "\tfeature\tmissingFeatureEntryFile1_" + transcriptComparison.getQueryTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                }
                                if (regionComparison.isMissingInQueryFile()) {
                                    writer.write(targetGeneId + "\t" + queryGeneId+ "\tfeature\tmissingFeatureEntryFile2_" + transcriptComparison.getQueryTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                }
                            }
                        }
                    }
                }
            }

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
