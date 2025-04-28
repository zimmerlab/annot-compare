package com.github.zimmerlab.gtfcompare.runner;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.zimmerlab.gtfcompare.GenomeSequenceExtractor;
import com.github.zimmerlab.gtfcompare.model.comparison.FeatureComparisonResult;
import com.github.zimmerlab.gtfcompare.model.comparison.RegionComparison;
import com.github.zimmerlab.gtfcompare.parser.FidxParser;
import compare.GTFCompare;
import org.apache.commons.cli.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.boot.CommandLineRunner;
import org.springframework.context.annotation.Profile;
import org.springframework.stereotype.Service;
import org.springframework.util.StopWatch;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

@Profile("analysis")
@Service
public class AnalysisRunner implements CommandLineRunner {

    private final static Logger LOG = LogManager.getLogger(TestRunner.class);

    public AnalysisRunner() {

    }

    @Override
    public void run(String... args) throws Exception {
        Options o = new Options();
        o.addOption(Option.builder()
                .option("h")
                .longOpt("help")
                .desc("Print the help message")
                .build());
        o.addOption(Option.builder()
                .longOpt("gtf")
                .numberOfArgs(1)
                .required()
                .desc("Path to gtf file")
                .type(File.class)
                .build());

        o.addOption(Option.builder()
                .longOpt("fasta")
                .numberOfArgs(1)
                .required()
                .desc("Path to fasta file")
                .type(File.class)
                .build());

        o.addOption(Option.builder()
                .longOpt("fidx")
                .numberOfArgs(1)
                .required()
                .desc("Path to fasta index file")
                .type(File.class)
                .build());

        o.addOption(Option.builder()
                .longOpt("gtf2")
                .numberOfArgs(1)
                .required()
                .desc("Path to gtf file")
                .type(File.class)
                .build());

        o.addOption(Option.builder()
                .longOpt("fasta2")
                .numberOfArgs(1)
                .required()
                .desc("Path to fasta file")
                .type(File.class)
                .build());

        o.addOption(Option.builder()
                .longOpt("fidx2")
                .numberOfArgs(1)
                .required()
                .desc("Path to fasta index file")
                .type(File.class)
                .build());

        o.addOption(Option.builder()
                .longOpt("o")
                .numberOfArgs(1)
                .required()
                .desc("Path to output file")
                .type(File.class)
                .build());

        CommandLineParser parser = new DefaultParser();

        CommandLine cmd = null;

        try {
            cmd = parser.parse(o, args);
        } catch (ParseException e) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Test", o, true);
            System.exit(1);
        }

        LOG.info("Running test");

        if (!cmd.hasOption("gtf")) {
            LOG.error("No gtf file specified");
            System.exit(1);
        }

        if (!cmd.hasOption("fasta")) {
            LOG.error("No fasta file specified");
            System.exit(1);
        }

        if (!cmd.hasOption("fidx")) {
            LOG.error("No fidx file specified");
            System.exit(1);
        }

        var fidxEntries = FidxParser.parse(cmd.getOptionValue("fidx"));
        var fidx2Entries = FidxParser.parse(cmd.getOptionValue("fidx2"));

        GtfFile gtfFile = new GtfFile(new File(cmd.getOptionValue("gtf")));
        GtfFile gtfFile2 = new GtfFile(new File(cmd.getOptionValue("gtf2")));

        var genomeSequenceExtractor1 = new GenomeSequenceExtractor(new File(cmd.getOptionValue("fasta")), fidxEntries);
        var genomeSequenceExtractor2 = new GenomeSequenceExtractor(new File(cmd.getOptionValue("fasta2")), fidx2Entries);

        try {
            LOG.info("Starting to parse gtfs");
            var stopWatch = new StopWatch();
            stopWatch.start("gtf parsing");
            gtfFile.parseAllContigs();
            LOG.info("GTF 1 parsed");
            gtfFile2.parseAllContigs();
            LOG.info("GTF 2 parsed");
            stopWatch.stop();
        } catch (java.text.ParseException e) {
            LOG.error("Error while parsing contigs", e);
            System.exit(-1);
        }


        var geneIds = gtfFile.getAllGeneFeatureIds();
        try (var writer = new BufferedWriter(new FileWriter(cmd.getOptionValue("o")))) {
            writer.write("geneId\tisSameGene\tdifference\n");

            geneLoop:
            for (String geneId : geneIds) {
                var geneFeature1 = gtfFile.getGeneFeature(geneId);
                var geneFeature2 = gtfFile2.getGeneFeature(geneId);

                if(geneFeature2 == null){
                    writer.write(geneId + "\tfalse\tgeneMissingInFile2\n");
                    continue geneLoop;
                }
                var comparisonResult = GTFCompare.compare(geneId, geneId, geneFeature1, geneFeature2, genomeSequenceExtractor1, genomeSequenceExtractor2);

                if (comparisonResult.areSameGene()) {
                    writer.write(geneId + "\ttrue\n");
                } else {
                    writer.write(geneId + "\tfalse\t");
                    var geneComparison = comparisonResult.getGeneComparison();

                    if (geneComparison.isStartDifferent()) {
                        writer.write("gene_start\n");
                        continue geneLoop;
                    }
                    if (geneComparison.isStopDifferent()) {
                        writer.write( "gene_stop\n");
                        continue geneLoop;
                    }
                    if (geneComparison.isStrandDifferent()) {
                        writer.write( "strand\n");
                        continue geneLoop;
                    }
                    if (!geneComparison.getSequenceComparison().isSameSequence()) {
                        writer.write( "gene_sequence\n");
                        continue geneLoop;
                    }

                    var transcriptComparisonResult = comparisonResult.getTranscriptComparisons();

                    for(var transcriptComparison : transcriptComparisonResult){
                        if (transcriptComparison.isStartDifferent()) {
                            writer.write("transcriptStart_" + transcriptComparison.getTranscriptId() + "\n");
                            continue geneLoop;
                        }
                        if (transcriptComparison.isStopDifferent()) {
                            writer.write("transcriptStop_" + transcriptComparison.getTranscriptId() + "\n");
                            continue geneLoop;
                        }
                        if (transcriptComparison.isTranscriptMissingInGene1()) {
                            writer.write("transcriptMissingInFile1_" + transcriptComparison.getTranscriptId() + "\n");
                            continue geneLoop;
                        }
                        if (transcriptComparison.isTranscriptMissingInGene2()) {
                            writer.write("transcriptMissingInFile2_" + transcriptComparison.getTranscriptId() + "\n");
                            continue geneLoop;
                        }

                        for(FeatureComparisonResult featureComparison : transcriptComparison.getFeatureComparisons()){
                            if(featureComparison.isMissingInTranscript1()){
                                writer.write("featureMissingInTranscript1_" + transcriptComparison.getTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                continue geneLoop;
                            }
                            if(featureComparison.isMissingInTranscript2()){
                                writer.write("featureMissingInTranscript2_" + transcriptComparison.getTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                continue geneLoop;
                            }

                            for(RegionComparison regionComparison : featureComparison.getRegionComparisons()){
                                if(regionComparison.isLengthDifferenceFound()){
                                    writer.write("regionLengthDifferent_" + transcriptComparison.getTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                    continue geneLoop;
                                }
                                if(regionComparison.isPositionDifferenceFound()){
                                    writer.write("regionPositionDifferent_" + transcriptComparison.getTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                    continue geneLoop;
                                }
                                if(regionComparison.isSequenceDifferenceFound()){
                                    writer.write("regionSequenceDifferent_" + transcriptComparison.getTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                    continue geneLoop;
                                }
                            }
                        }
                    }
                }
            }

            for (String geneId : gtfFile2.getAllGeneFeatureIds()) {
                if(gtfFile.getGeneFeature(geneId) == null){
                    writer.write(geneId + "\tfalse\tgeneMissingInFile1\n");
                }
            }

        } catch (Exception e) {
            LOG.error("Error while comparing GTFs", e);
        }


    }
}