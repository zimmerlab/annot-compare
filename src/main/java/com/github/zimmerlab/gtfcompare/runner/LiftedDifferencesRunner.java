package com.github.zimmerlab.gtfcompare.runner;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.zimmerlab.gtfcompare.analysis.clique.CliqueAnalysisGenes;
import com.github.zimmerlab.gtfcompare.analysis.lifted.LiftedDifferences;
import com.github.zimmerlab.gtfcompare.parser.FidxParser;
import com.github.zimmerlab.gtfcompare.utils.GenomeSequenceExtractor;
import org.apache.commons.cli.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.boot.CommandLineRunner;
import org.springframework.context.annotation.Profile;
import org.springframework.stereotype.Service;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Objects;

@Profile("liftedDifferences")
@Service
public class LiftedDifferencesRunner implements CommandLineRunner {
    private final static Logger logger = LogManager.getLogger(LiftedDifferencesRunner.class);

    @Override
    public void run(String... args) throws Exception {
        Options o = new Options();
        o.addOption(Option.builder().option("h").longOpt("help").desc("Print the help message").build());
        o.addOption(Option.builder().longOpt("target-gtf").numberOfArgs(1).required().desc("Path to target gtf file").type(File.class).build());
        o.addOption(Option.builder().longOpt("query-gtf").numberOfArgs(1).required().desc("Path to query gtf file").type(File.class).build());
        o.addOption(Option.builder().longOpt("lifted-query-gtf").numberOfArgs(1).required().desc("Path to lifted query gtf file").type(File.class).build());
        o.addOption(Option.builder().longOpt("query-fasta").numberOfArgs(1).required().desc("Path to fasta file").type(File.class).build());
        o.addOption(Option.builder().longOpt("query-fidx").numberOfArgs(1).required().desc("Path to fasta index file").type(File.class).build());
        o.addOption(Option.builder().longOpt("target-fasta").numberOfArgs(1).required().desc("Path to fasta file").type(File.class).build());
        o.addOption(Option.builder().longOpt("target-fidx").numberOfArgs(1).required().desc("Path to fasta index file").type(File.class).build());
        o.addOption(Option.builder().longOpt("output").numberOfArgs(1).required().desc("Path to output file").type(File.class).build());

        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = null;

        try {
            cmd = parser.parse(o, args);
        } catch (ParseException e) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("benchmark", o, true);
            System.exit(1);
        }

        if (!cmd.hasOption("target-gtf")) {
            logger.error("No target gtf specified");
            System.exit(1);
        }

        if (!cmd.hasOption("query-gtf")) {
            logger.error("No query gtf specified");
            System.exit(1);
        }

        if (!cmd.hasOption("lifted-query-gtf")) {
            logger.error("No lifted query gtf specified");
            System.exit(1);
        }

        if (!cmd.hasOption("output")) {
            logger.error("No output specified");
            System.exit(1);
        }

        if (!cmd.hasOption("target-fasta")) {
            logger.error("No target fasta specified");
            System.exit(1);
        }

        if (!cmd.hasOption("target-fidx")) {
            logger.error("No target fidx specified");
            System.exit(1);
        }

        if (!cmd.hasOption("query-fasta")) {
            logger.error("No query fasta specified");
            System.exit(1);
        }

        if (!cmd.hasOption("query-fidx")) {
            logger.error("No query fidx specified");
            System.exit(1);
        }

        var queryPath = cmd.getOptionValue("query-gtf");
        var liftedQueryPath = cmd.getOptionValue("lifted-query-gtf");
        var targetPath = cmd.getOptionValue("target-gtf");

        var targetFasta = cmd.getOptionValue("target-fasta");
        var targetFidx = cmd.getOptionValue("target-fidx");

        var queryFasta = cmd.getOptionValue("query-fasta");
        var queryFidx = cmd.getOptionValue("query-fidx");

        var outputPath = cmd.getOptionValue("output");

        var targetGtf = new GtfFile(new File(targetPath));
        var queryGtf = new GtfFile(new File(queryPath));
        var liftedQueryGtf = new GtfFile(new File(liftedQueryPath));

        var queryFidxMapping = FidxParser.parse(queryFidx);
        logger.info("Query FIDX Mapping Size: {}", queryFidxMapping.size());
        logger.info("Creating Query GSE");
        var queryGSE = new GenomeSequenceExtractor(new File(queryFasta), queryFidxMapping);

        var targetFidxMapping = FidxParser.parse(targetFidx);
        logger.info("Target FIDX Mapping Size: {}", targetFidxMapping.size());
        logger.info("Creating Target GSE");
        var targetGSE = new GenomeSequenceExtractor(new File(targetFasta),  targetFidxMapping);
        try {
            logger.info("Starting analysis");
            while (true) {
                targetGtf.parseNextContig();
                queryGtf.parseNextContig();
                liftedQueryGtf.parseNextContig();

                String t = targetGtf.getParsedContig();
                String q = queryGtf.getParsedContig();
                String lq = liftedQueryGtf.getParsedContig();
                logger.info("Using target: {}, query: {}, lifted query: {}",t, q, lq);
                if (!Objects.equals(t, q) || !Objects.equals(q, lq))
                    throw new Exception("Contigs do not match Target: " + t + ", Query: " + q + ", Lifted Query: " + lq);

                LiftedDifferences.analyze(t, queryGtf, targetGtf, liftedQueryGtf, queryGSE, targetGSE, outputPath);
            }
        } catch (java.text.ParseException e) {
            logger.info("Program finished: {}", e.getMessage());
        } catch (Exception e) {
            logger.error("Program failed", e);
        }

        LiftedDifferences.printResults(outputPath);

    }
}
