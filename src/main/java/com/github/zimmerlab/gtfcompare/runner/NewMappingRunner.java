package com.github.zimmerlab.gtfcompare.runner;


import com.github.kleinsamuel.gtfutils.GtfConfig;
import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.zimmerlab.gtfcompare.newmapping.Mapping;
import com.github.zimmerlab.gtfcompare.newmapping.outpututil.MappingOutputWriter;
import com.github.zimmerlab.gtfcompare.newmapping.outpututil.UnmappedWriter;
import com.github.zimmerlab.gtfcompare.parser.FidxParser;
import com.github.zimmerlab.gtfcompare.utils.GenomeSequenceExtractor;
import org.apache.commons.cli.*;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.boot.CommandLineRunner;
import org.springframework.context.annotation.Profile;
import org.springframework.stereotype.Service;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.util.HashSet;
import java.util.Objects;
import java.util.Set;

@Profile("newMapping")
@Service
public class NewMappingRunner implements CommandLineRunner {
    private final static Logger logger = LogManager.getLogger(NewMappingRunner.class);

    @Override
    public void run(String... args) throws Exception {
        Options o = new Options();
        o.addOption(Option.builder().option("h").longOpt("help").desc("Print the help message").build());
        o.addOption(Option.builder().longOpt("target-gtf").numberOfArgs(1).required().desc("Path to target gtf file").type(File.class).build());
        o.addOption(Option.builder().longOpt("target-fasta").numberOfArgs(1).required().desc("Path to target fasta file").type(File.class).build());
        o.addOption(Option.builder().longOpt("query-fasta").numberOfArgs(1).required().desc("Path to query fasta file").type(File.class).build());
        o.addOption(Option.builder().longOpt("query-fai").numberOfArgs(1).required().desc("Path to query fai file").type(File.class).build());
        o.addOption(Option.builder().longOpt("target-fai").numberOfArgs(1).required().desc("Path to target fai file").type(File.class).build());
        o.addOption(Option.builder().longOpt("query-gtf").numberOfArgs(1).required().desc("Path to query gtf file").type(File.class).build());
        o.addOption(Option.builder().longOpt("output").numberOfArgs(1).required().desc("Path to output file").type(File.class).build());
        o.addOption(Option.builder().longOpt("allowed-types").hasArg().desc("Comma-separated list of allowed types").build());
        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = null;

        try {
            cmd = parser.parse(o, args);
        } catch (ParseException e) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("newMapping", o, true);
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

        if (!cmd.hasOption("output")) {
            logger.error("No output specified");
            System.exit(1);
        }

        Set<String> allowedTypes = null;

        if (cmd.hasOption("allowed-types")) {
            var set = new HashSet<String>();

            for (var type : cmd.getOptionValue("allowed-types").split(",")) {
                var defaultType = GtfConfig.getDefault(type);
                if (defaultType != null) {
                    set.add(defaultType);
                } else {
                    logger.info("Skipping type {} - no default value found", type);
                }
            }

            allowedTypes = set.isEmpty() ? null : set;
        }

        if (cmd.hasOption("allowed-types") && allowedTypes == null) {
            logger.info("no valid types given - using default values for allowed-types");
        }


        var queryPath = cmd.getOptionValue("query-gtf");
        var targetPath = cmd.getOptionValue("target-gtf");
        var targetFastaPath = new File(cmd.getOptionValue("target-fasta"));
        var queryFastaPath = new File(cmd.getOptionValue("query-fasta"));
        var targetFaiPath = cmd.getOptionValue("target-fai");
        var queryFaiPath = cmd.getOptionValue("query-fai");
        var outputPath = cmd.getOptionValue("output");

        var targetGtf = new GtfFile(new File(targetPath));
        var queryGtf = new GtfFile(new File(queryPath));

        var targetFai = FidxParser.parse(targetFaiPath);
        var queryFai = FidxParser.parse(queryFaiPath);

        var targetSequenceExtractor = new GenomeSequenceExtractor(targetFastaPath, targetFai);
        var querySequenceExtractor = new GenomeSequenceExtractor(queryFastaPath, queryFai);


        try (var writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputPath), StandardCharsets.UTF_8));
             var unmappedWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputPath + ".unmapped"), StandardCharsets.UTF_8))) {

            writer.write("contig\tqueryId\ttargetId\tqueryTranscriptId\ttargetTranscriptId\tmapping_origins\n");
            unmappedWriter.write("contig\tgeneId\tname\torigin\n");
            while (true) {
                targetGtf.parseNextContig();
                queryGtf.parseNextContig();

                String t = targetGtf.getParsedContig();
                String q = queryGtf.getParsedContig();
                if (!Objects.equals(t, q)) throw new Exception("Contigs do not match");

                logger.info("Current Contig: {}", t);
                var res = Mapping.map(targetGtf, queryGtf, allowedTypes, targetSequenceExtractor, querySequenceExtractor);
                MappingOutputWriter.write(res.results(), t, writer);
                UnmappedWriter.write(res.unmappedQueries(), q, "QUERY", unmappedWriter);
                UnmappedWriter.write(res.unmappedTargets(), t, "TARGET", unmappedWriter);
            }
        } catch (java.text.ParseException e) {
            logger.debug("Program finished: {}", e.getMessage());
        } catch (Exception e) {
            logger.error("Program failed", e);
        }

    }
}

