package com.github.zimmerlab.gtfcompare.runner;


import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.zimmerlab.gtfcompare.analysis.clique.CliqueAnalysisGenes;
import com.github.zimmerlab.gtfcompare.analysis.clique.Reporter;
import com.github.zimmerlab.gtfcompare.newmapping.Mapping;
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

@Profile("newMapping")
@Service
public class NewMappingRunner implements CommandLineRunner {
    private final static Logger logger = LogManager.getLogger(NewMappingRunner.class);
    @Override
    public void run(String... args) throws Exception {
        Options o = new Options();
        o.addOption(Option.builder().option("h").longOpt("help").desc("Print the help message").build());
        o.addOption(Option.builder().longOpt("target-gtf").numberOfArgs(1).required().desc("Path to target gtf file").type(File.class).build());
        o.addOption(Option.builder().longOpt("query-gtf").numberOfArgs(1).required().desc("Path to query gtf file").type(File.class).build());
        o.addOption(Option.builder().longOpt("output").numberOfArgs(1).required().desc("Path to output file").type(File.class).build());

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

        var queryPath = cmd.getOptionValue("query-gtf");
        var targetPath = cmd.getOptionValue("target-gtf");
        var outputPath = cmd.getOptionValue("output");

        var targetGtf = new GtfFile(new File(targetPath));
        var queryGtf = new GtfFile(new File(queryPath));


        try {
            while (true) {
                targetGtf.parseNextContig();
                queryGtf.parseNextContig();

                String t = targetGtf.getParsedContig();
                String q = queryGtf.getParsedContig();
                if (!Objects.equals(t, q)) throw new Exception("Contigs do not match");

                var res = Mapping.map(targetGtf, queryGtf);
            }
        } catch (java.text.ParseException e) {
            logger.debug("Program finished: {}", e.getMessage());
        } catch (Exception e) {
            logger.error("Program failed", e);
        }

    }
}
