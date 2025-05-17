package com.github.zimmerlab.gtfcompare.runner;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.zimmerlab.gtfcompare.AnnotComparator;
import com.github.zimmerlab.gtfcompare.compare.ComparisonConfig;
import com.github.zimmerlab.gtfcompare.compare.ComparisonConfigBuilder;
import com.github.zimmerlab.gtfcompare.model.config.ConfigJSON;
import com.github.zimmerlab.gtfcompare.parser.FidxParser;
import com.github.zimmerlab.gtfcompare.utils.Constants;
import com.github.zimmerlab.gtfcompare.utils.GenomeSequenceExtractor;
import org.apache.commons.cli.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.boot.CommandLineRunner;
import org.springframework.context.annotation.Profile;
import org.springframework.stereotype.Service;
import org.springframework.util.StopWatch;

import java.io.File;
import java.util.*;

@Profile("analysis")
@Service
public class AnalysisRunner implements CommandLineRunner {
    private final static Logger LOG = LogManager.getLogger(AnalysisRunner.class);
    private final static List<StopWatch> stopWatches = Constants.STOP_WATCHES;

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

        o.addOption(Option.builder()
                .longOpt("config")
                .numberOfArgs(1)
                .required()
                .desc("Path to config file")
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

        if (!cmd.hasOption("config")) {
            LOG.error("No config file specified");
            System.exit(1);
        }

        if (!cmd.hasOption("o")) {
            LOG.error("No output path specified");
            System.exit(1);
        }

        var fidxEntries = FidxParser.parse(cmd.getOptionValue("fidx"));
        var fidx2Entries = FidxParser.parse(cmd.getOptionValue("fidx2"));

        GtfFile gtfFile = new GtfFile(new File(cmd.getOptionValue("gtf")));
        GtfFile gtfFile2 = new GtfFile(new File(cmd.getOptionValue("gtf2")));

        var targetSequenceExtractor = new GenomeSequenceExtractor(new File(cmd.getOptionValue("fasta")), fidxEntries);
        var querySequenceExtractor = new GenomeSequenceExtractor(new File(cmd.getOptionValue("fasta2")), fidx2Entries);

        gtfFile.parseAllContigs();
        gtfFile2.parseAllContigs();

        var configPath = cmd.getOptionValue("config");

        var objectMapper = new ObjectMapper();
        ConfigJSON jsonConfig = null;
        try {
            jsonConfig = objectMapper.readValue(new File(configPath), ConfigJSON.class);
        } catch (Exception e) {
            LOG.error("Error reading config file: {}", e.getMessage());
            System.exit(1);
        }

        var config = getComparisonConfig(jsonConfig);

        var annotCompare = new AnnotComparator(gtfFile, gtfFile2, targetSequenceExtractor, querySequenceExtractor, config, cmd.getOptionValue("o"));
        annotCompare.compare();
    }
    private static ComparisonConfig getComparisonConfig(ConfigJSON jsonConfig) {
        var configBuilder = new ComparisonConfigBuilder();

        var configFeatures = jsonConfig.getFeatures();

        var length = configFeatures.get(Constants.LENGTH_COMPARATOR_NAME);
        if (length != null && length.isEnabled()) {
            configBuilder.enableFeature(Constants.LENGTH_COMPARATOR_NAME);
            var th = length.getThreshold();
            if (th != null) {
                configBuilder.setThreshold(Constants.LENGTH_COMPARATOR_NAME, th);
            }
        }

        var start = configFeatures.get(Constants.START_COMPARATOR_NAME);
        if (start != null && start.isEnabled()) {
            configBuilder.enableFeature(Constants.START_COMPARATOR_NAME);
            var th = start.getThreshold();
            if (th != null) {
                configBuilder.setThreshold(Constants.START_COMPARATOR_NAME, th);
            }
        }

        var stop = configFeatures.get(Constants.STOP_COMPARATOR_NAME);
        if (stop != null && stop.isEnabled()) {
            configBuilder.enableFeature(Constants.STOP_COMPARATOR_NAME);
            var th = stop.getThreshold();
            if (th != null) {
                configBuilder.setThreshold(Constants.STOP_COMPARATOR_NAME, th);
            }
        }

        var seq = configFeatures.get(Constants.SEQUENCE_COMPARATOR_NAME);
        if (seq != null && seq.isEnabled()) {
            configBuilder.enableFeature(Constants.SEQUENCE_COMPARATOR_NAME);
        }

        return configBuilder.build();
    }
}
