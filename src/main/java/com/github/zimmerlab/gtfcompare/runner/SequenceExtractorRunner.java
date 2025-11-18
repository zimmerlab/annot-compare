package com.github.zimmerlab.gtfcompare.runner;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.zimmerlab.gtfcompare.parser.FidxParser;
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


        if (!cmd.hasOption("o")) {
            LOG.error("No output path specified");
            System.exit(1);
        }




        var fidxEntries = FidxParser.parse(cmd.getOptionValue("fidx"));

        GtfFile gtfFile = new GtfFile(new File(cmd.getOptionValue("gtf")));

        var targetSequenceExtractor = new GenomeSequenceExtractor(new File(cmd.getOptionValue("fasta")), fidxEntries);

        gtfFile.parseAllContigs();

        for(var geneId : gtfFile.getAllGeneFeatureIds()){
            for(var transcript : gtfFile.getGeneFeature(geneId).getTranscripts()){
                var transcriptSeq = targetSequenceExtractor.getSequence(transcript.getBaseData().getContig(), transcript.getBaseData().getStart(), transcript.getBaseData().getEnd());
                try (BufferedWriter writer = Files.newBufferedWriter(Path.of(cmd.getOptionValue("o"), transcript.getTranscriptId() + ".fasta"))) {
                    writer.write(">" + transcript.getTranscriptId() + "\n");
                    writer.write(transcriptSeq + "\n");

                }
            }
        }



    }

}
