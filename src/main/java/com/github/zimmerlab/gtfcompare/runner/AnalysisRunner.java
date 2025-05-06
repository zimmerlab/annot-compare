package com.github.zimmerlab.gtfcompare.runner;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.zimmerlab.gtfcompare.utils.GenomeSequenceExtractor;
import com.github.zimmerlab.gtfcompare.model.comparison.FeatureComparisonResult;
import com.github.zimmerlab.gtfcompare.model.comparison.RegionComparison;
import com.github.zimmerlab.gtfcompare.parser.FidxParser;
import com.github.zimmerlab.gtfcompare.utils.Constants;
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
import java.util.List;

@Profile("analysis")
@Service
public class AnalysisRunner implements CommandLineRunner {

    private final static Logger LOG = LogManager.getLogger(TestRunner.class);
    private final static List<StopWatch> stopWatches = Constants.STOP_WATCHES;
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
            stopWatches.add(stopWatch);
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
            writer.write("geneId\tcategory\tdifference\n");

            for (String geneId : geneIds) {
                var geneFeature1 = gtfFile.getGeneFeature(geneId);
                var geneFeature2 = gtfFile2.getGeneFeature(geneId);

                if(geneFeature2 == null){
                    writer.write(geneId + "\tgene\tmissingInFile2\n");
                    continue;
                }

                var comparisonResult = GTFCompare.compare(geneId, geneId, geneFeature1, geneFeature2, genomeSequenceExtractor1, genomeSequenceExtractor2);

                if (comparisonResult.areSameGene()) {
                    writer.write(geneId + "\n");
                } else {
                    var geneComparison = comparisonResult.getGeneComparison();

                    if (geneComparison.isStartDifferent()) {
                        writer.write(geneId + "\tgene\tstart\n");
                    }
                    if (geneComparison.isStopDifferent()) {
                        writer.write( geneId + "\tgene\tstop\n");
                    }
                    if (geneComparison.isStrandDifferent()) {
                        writer.write(geneId + "\tgene\tstrand\n");
                    }
                    if (!geneComparison.getSequenceComparison().isSameSequence()) {
                        writer.write( geneId + "\tgene\tseq\n");
                    }
                    if(geneComparison.isDifferentLength()){
                        writer.write( geneId + "\tgene\tlength\n");
                    }

                    var transcriptComparisonResult = comparisonResult.getTranscriptComparisons();

                    for(var transcriptComparison : transcriptComparisonResult){
                        if (transcriptComparison.isStartDifferent()) {
                            writer.write(geneId + "\ttranscript\tstart_" + transcriptComparison.getTranscriptId() + "\n");
                        }
                        if (transcriptComparison.isStopDifferent()) {
                            writer.write(geneId + "\ttranscript\tstop_" + transcriptComparison.getTranscriptId() + "\n");
                        }
                        if(transcriptComparison.isSequenceDifferent()){
                            writer.write(geneId + "\ttranscript\tseq_" + transcriptComparison.getTranscriptId() + "\n");
                        }
                        if (transcriptComparison.isTranscriptMissingInGene1()) {
                            writer.write(geneId + "\ttranscript\tmissingInFile1_" + transcriptComparison.getTranscriptId() + "\n");
                        }
                        if (transcriptComparison.isTranscriptMissingInGene2()) {
                            writer.write(geneId + "\ttranscript\tmissingInFile2_" + transcriptComparison.getTranscriptId() + "\n");
                        }
                        if(transcriptComparison.isLengthDifferent()){
                            writer.write(geneId + "\ttranscript\tlength_" + transcriptComparison.getTranscriptId() + "\n");

                        }

                        for(FeatureComparisonResult featureComparison : transcriptComparison.getFeatureComparisons()){
                            if(featureComparison.isMissingInTranscript1()){
                                writer.write(geneId + "\tfeature\tmissingInTranscript1_" + transcriptComparison.getTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                            }
                            if(featureComparison.isMissingInTranscript2()){
                                writer.write(geneId + "\tfeature\tmissingInTranscript2_" + transcriptComparison.getTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                            }

                            for(RegionComparison regionComparison : featureComparison.getRegionComparisons()){
                                if(regionComparison.isLengthDifferenceFound()){
                                    writer.write(geneId + "\tfeature\tlength_" + transcriptComparison.getTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                }
                                if(regionComparison.isPositionDifferenceFound()){
                                    writer.write(geneId + "\tfeature\tposition_" + transcriptComparison.getTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                }
                                if(regionComparison.isSequenceDifferenceFound()){
                                    writer.write(geneId + "\tfeature\tseq_" + transcriptComparison.getTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                }
                                if(regionComparison.isMissingInFile1()){
                                    writer.write(geneId + "\tfeature\tmissingFeatureEntryFile1_" + transcriptComparison.getTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                }
                                if(regionComparison.isMissingInFile2()){
                                    writer.write(geneId + "\tfeature\tmissingFeatureEntryFile2_" + transcriptComparison.getTranscriptId() + "_" + featureComparison.getFeatureType() + "\n");
                                }
                            }
                        }
                    }
                }
            }

            for (String geneId : gtfFile2.getAllGeneFeatureIds()) {
                if(gtfFile.getGeneFeature(geneId) == null){
                    writer.write(geneId + "\tgene\tmissingInFile1\n");
                }
            }

        } catch (Exception e) {
            LOG.error("Error while comparing GTFs", e);
        }


    }
}