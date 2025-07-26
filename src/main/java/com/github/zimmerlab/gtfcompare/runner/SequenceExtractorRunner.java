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

import java.io.BufferedWriter;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;

@Profile("seqExtractor")
@Service
public class SequenceExtractorRunner implements CommandLineRunner {
    private final static Logger LOG = LogManager.getLogger(SequenceExtractorRunner.class);

    public SequenceExtractorRunner() {
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

        if (!cmd.hasOption("fasta2")) {
            LOG.error("No fasta2 file specified");
            System.exit(1);
        }

        if (!cmd.hasOption("fidx")) {
            LOG.error("No fidx file specified");
            System.exit(1);
        }

        if (!cmd.hasOption("fidx2")) {
            LOG.error("No fidx2 file specified");
            System.exit(1);
        }

        if (!cmd.hasOption("o")) {
            LOG.error("No output path specified");
            System.exit(1);
        }

        var gene1 = "ENSG00000199038";
        var gene2 = "ENSG00000254474";

        var fidxEntries = FidxParser.parse(cmd.getOptionValue("fidx"));
        var fidx2Entries = FidxParser.parse(cmd.getOptionValue("fidx2"));

        GtfFile gtfFile = new GtfFile(new File(cmd.getOptionValue("gtf")));
        GtfFile gtfFile2 = new GtfFile(new File(cmd.getOptionValue("gtf2")));

        var targetSequenceExtractor = new GenomeSequenceExtractor(new File(cmd.getOptionValue("fasta")), fidxEntries);
        var querySequenceExtractor = new GenomeSequenceExtractor(new File(cmd.getOptionValue("fasta2")), fidx2Entries);

        gtfFile.parseAllContigs();
        gtfFile2.parseAllContigs();

        var targetGene1 = gtfFile.getGeneFeature(gene1);
        var queryGene1 = gtfFile2.getGeneFeature(gene1);

        var targetGene2 = gtfFile.getGeneFeature(gene2);
        var queryGene2 = gtfFile2.getGeneFeature(gene2);

        var targetGene1Seq = targetSequenceExtractor.getSequence(targetGene1.getBaseData().getContig(), targetGene1.getBaseData().getStart(), targetGene1.getBaseData().getEnd());
        var queryGene1Seq = querySequenceExtractor.getSequence(queryGene1.getBaseData().getContig(), queryGene1.getBaseData().getStart(), queryGene1.getBaseData().getEnd());

        var targetGene2Seq = targetSequenceExtractor.getSequence(targetGene2.getBaseData().getContig(), targetGene2.getBaseData().getStart(), targetGene2.getBaseData().getEnd());
        var queryGene2Seq = querySequenceExtractor.getSequence(queryGene2.getBaseData().getContig(), queryGene2.getBaseData().getStart(), queryGene2.getBaseData().getEnd());

        System.out.println(targetGene1Seq);
        System.out.println(targetGene2Seq);
        System.out.println(queryGene1Seq);
        System.out.println(queryGene2Seq);

        try (BufferedWriter writer = Files.newBufferedWriter(Path.of(cmd.getOptionValue("o"), "71_72_" + gene1 + "_" + gene2 + ".fasta"))) {
            writer.write(">" + targetGene1.getGeneId() + "|target1\n");
            writer.write(targetGene1Seq + "\n");
            writer.write(">" + queryGene1.getGeneId() + "|query1\n");
            writer.write(queryGene1Seq + "\n");
            writer.write(">" + targetGene2.getGeneId() + "|target2\n");
            writer.write(targetGene2Seq + "\n");
            writer.write(">" + queryGene2.getGeneId() + "|query2\n");
            writer.write(queryGene2Seq + "\n");

        }
    }

}
