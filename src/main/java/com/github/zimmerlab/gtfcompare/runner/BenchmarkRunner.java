package com.github.zimmerlab.gtfcompare.runner;

import com.github.zimmerlab.gtfcompare.utils.JfrAggregate;
import org.apache.commons.cli.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.boot.CommandLineRunner;
import org.springframework.context.annotation.Profile;
import org.springframework.stereotype.Service;

import java.io.File;

@Profile("benchmark")
@Service
public class BenchmarkRunner implements CommandLineRunner {
    private final static Logger logger = LogManager.getLogger(BenchmarkRunner.class);
    @Override
    public void run(String... args) throws Exception{
        Options o = new Options();
        o.addOption(Option.builder().option("h").longOpt("help").desc("Print the help message").build());
        o.addOption(Option.builder().longOpt("jfr").numberOfArgs(1).required().desc("Path to jfr file").type(File.class).build());
        o.addOption(Option.builder().longOpt("o").numberOfArgs(1).required().desc("Path to output directory").type(File.class).build());
        o.addOption(Option.builder().longOpt("pre").numberOfArgs(1).required().desc("Prefix for output files").type(File.class).build());

        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = null;

        try {
            cmd = parser.parse(o, args);
        } catch (ParseException e) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Test", o, true);
            System.exit(1);
        }

        if (!cmd.hasOption("jfr")) {
            logger.error("No gtf file specified");
            System.exit(1);
        }

        if (!cmd.hasOption("o")) {
            logger.error("No output directory specified");
            System.exit(1);
        }

        if (!cmd.hasOption("pre")) {
            logger.error("No prefix for output files");
            System.exit(1);
        }

        JfrAggregate.benchmark(cmd.getOptionValue("jfr"), cmd.getOptionValue("o"), cmd.getOptionValue("pre"));
    }


}
