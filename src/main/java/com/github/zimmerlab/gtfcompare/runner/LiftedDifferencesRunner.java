package com.github.zimmerlab.gtfcompare.runner;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.zimmerlab.gtfcompare.analysis.lifted.LiftedDifferences;
import org.apache.commons.cli.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.boot.CommandLineRunner;
import org.springframework.context.annotation.Profile;
import org.springframework.stereotype.Service;

import java.io.File;

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

        var queryPath = cmd.getOptionValue("query-gtf");
        var liftedQueryPath = cmd.getOptionValue("lifted-query-gtf");
        var targetPath = cmd.getOptionValue("target-gtf");
        var outputPath = cmd.getOptionValue("output");

        var targetGtf = new GtfFile(new File(targetPath));
        var queryGtf = new GtfFile(new File(queryPath));
        var liftedQueryGtf = new GtfFile(new File(liftedQueryPath));

        LiftedDifferences.analyze(queryGtf, targetGtf, liftedQueryGtf);
    }
}
