package com.github.zimmerlab.gtfcompare;

import com.github.kleinsamuel.gtfutils.GtfConfig;
import com.github.kleinsamuel.gtfutils.GtfConstants;
import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.GeneFeature;
import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.parser.FidxParser;
import com.github.zimmerlab.gtfcompare.*;
import compare.GTFCompare;
import org.apache.commons.cli.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.boot.CommandLineRunner;
import org.springframework.context.annotation.Profile;
import org.springframework.stereotype.Service;

import java.io.File;
import java.util.List;
import java.util.Set;

@Profile("test")
@Service
public class TestRunner implements CommandLineRunner {

    private final static Logger LOG = LogManager.getLogger(TestRunner.class);

    public TestRunner() {

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

        var genomeSequenceExtractor = new GenomeSequenceExtractor(new File(cmd.getOptionValue("fasta")), fidxEntries);
        var genomeSequenceExtractor2 = new GenomeSequenceExtractor(new File(cmd.getOptionValue("fasta2")), fidxEntries);
        GtfFile gtfFile = new GtfFile(new File(cmd.getOptionValue("gtf")));
        GtfFile gtfFile2 = new GtfFile(new File(cmd.getOptionValue("gtf2")));

        // parse all contigs from the gtf file
//        gtfFile.parseAllContigs();

        // parse the next (first in this case) contig (chr) from the gtf file
        do {
            gtfFile.parseNextContig();
        }  while(!gtfFile.getParsedContig().equals("3"));

        LOG.info("Parsed GTF 1");
        do{
            gtfFile2.parseNextContig();
        } while (!gtfFile2.getParsedContig().equals("3"));
        LOG.info("Parsed GTF 2");

        LOG.info("parsed contig: {}", gtfFile.getParsedContig());

        var geneIds = gtfFile.getAllGeneFeatureIds();

        for(String geneId : geneIds){
            var geneFeature1 = gtfFile.getGeneFeature(geneId);
            var geneFeature2 = gtfFile2.getGeneFeature(geneId);
            GTFCompare.comparePosition(geneFeature1, geneFeature2);
        }

    }
}
