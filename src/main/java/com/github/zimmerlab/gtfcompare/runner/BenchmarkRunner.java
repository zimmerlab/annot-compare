package com.github.zimmerlab.gtfcompare.runner;

import com.github.zimmerlab.gtfcompare.utils.JfrAggregate;
import org.apache.commons.cli.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.boot.CommandLineRunner;
import org.springframework.context.annotation.Profile;
import org.springframework.stereotype.Service;

import java.io.File;
import java.nio.file.Files;

@Profile("benchmark")
@Service
public class BenchmarkRunner implements CommandLineRunner {
    private final static Logger logger = LogManager.getLogger(BenchmarkRunner.class);
    @Override
    public void run(String... args) throws Exception{
        Options o = new Options();
        o.addOption(Option.builder().option("h").longOpt("help").desc("Print the help message").build());
        o.addOption(Option.builder().longOpt("dir").numberOfArgs(1).required().desc("Path to jfr folder").type(File.class).build());
        o.addOption(Option.builder().longOpt("o").numberOfArgs(1).required().desc("Path to output directory").type(File.class).build());
        o.addOption(Option.builder().longOpt("pre").numberOfArgs(1).desc("Prefix for output files").type(File.class).build());

        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = null;

        try {
            cmd = parser.parse(o, args);
        } catch (ParseException e) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("benchmark", o, true);
            System.exit(1);
        }

        if (!cmd.hasOption("dir")) {
            logger.error("No jfr directory specified");
            System.exit(1);
        }

        if (!cmd.hasOption("o")) {
            logger.error("No output directory specified");
            System.exit(1);
        }

        String rootDir = cmd.getOptionValue("dir");
        String outDir = cmd.getOptionValue("o");
        String prefix = cmd.getOptionValue("pre");

        Files.createDirectories(java.nio.file.Paths.get(outDir));

        var root = java.nio.file.Paths.get(rootDir);
        if (!java.nio.file.Files.isDirectory(root)) {
            logger.error("--dir {} is not a directory", rootDir);
            System.exit(1);
        }

        try (var subs = Files.list(root)) {
            subs.filter(Files::isDirectory).forEach(sub -> {
                try {
                    try (var files = java.nio.file.Files.list(sub)) {
                        var jfrOpt = files
                                .filter(f -> f.getFileName().toString().endsWith(".jfr"))
                                .findFirst();
                        if (jfrOpt.isEmpty()) {
                            logger.warn("No JFR file found in subdir {}", sub);
                            return;
                        }
                        var jfrFile = jfrOpt.get().toString();
                        var subPrefix = prefix != null ?  prefix + "_" + sub.getFileName().toString() : sub.getFileName().toString();
                        JfrAggregate.benchmark(jfrFile, outDir, subPrefix);
                        logger.info("Processed JFR {} with prefix {}", jfrFile, subPrefix);
                    }
                } catch (Exception ex) {
                    logger.error("Failed to process subdir {}: {}", sub, ex.getMessage(), ex);
                }
            });
        }
    }


}
