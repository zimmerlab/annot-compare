package com.github.zimmerlab.gtfcompare.runner;

import com.github.zimmerlab.gtfcompare.newmapping.FilterMappings;
import com.github.zimmerlab.gtfcompare.newmapping.model.MappingOrigin;
import org.apache.commons.cli.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.boot.CommandLineRunner;
import org.springframework.context.annotation.Profile;
import org.springframework.stereotype.Service;

import java.io.File;
import java.util.HashSet;

@Profile("filterGeneMappings")
@Service
public class FilterMappingsRunner implements CommandLineRunner {
    private final static Logger logger = LogManager.getLogger(FilterMappingsRunner.class);
    @Override
    public void run(String... args) throws Exception {
        Options o = new Options();
        o.addOption(Option.builder().longOpt("mapping-file").numberOfArgs(1).required().desc("Path to mapping file").type(File.class).build());
        o.addOption(Option.builder().longOpt("output-file").numberOfArgs(1).required().desc("Path to output file").type(File.class).build());
        o.addOption(Option.builder().longOpt("require-all").numberOfArgs(1).desc("Require all allowed mappings").type(Boolean.class).build());
        o.addOption(Option.builder().longOpt("allowed-mappings").hasArg().desc("Comma-separated list of allowed mappings").build());
        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = null;

        try {
            cmd = parser.parse(o, args);
        } catch (ParseException e) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("filterGeneMappings", o, true);
            System.exit(1);
        }

        if (!cmd.hasOption("mapping-file")) {
            logger.error("No mapping file specified");
            System.exit(1);
        }

        if (!cmd.hasOption("output-file")) {
            logger.error("No output file specified");
            System.exit(1);
        }

        if (!cmd.hasOption("allowed-mappings")) {
            logger.error("No allowed mappings specified");
            System.exit(1);
        }

        boolean requireAll = true;

        if (cmd.hasOption("require-all")) {
            requireAll = Boolean.parseBoolean(cmd.getOptionValue("require-all"));
        }

        var mappingFile = cmd.getOptionValue("mapping-file");
        var allowedMappingsString = cmd.getOptionValue("allowed-mappings").split(",");
        var allowedMappings = new HashSet<MappingOrigin>();
        for(var allowedMapping : allowedMappingsString){
            switch (allowedMapping.trim().toLowerCase()){
                case "gene-id":
                    allowedMappings.add(MappingOrigin.GENE_ID_MAPPING);
                    break;
                case "gene-name":
                    allowedMappings.add(MappingOrigin.NAME_MAPPING);
                    break;
                case "transcript-id":
                    allowedMappings.add(MappingOrigin.TRANSCRIPT_ID_MAPPING);
                    break;
                case "overlapping":
                    allowedMappings.add(MappingOrigin.OVERLAPPING);
                    break;
                case "distance":
                    allowedMappings.add(MappingOrigin.DISTANCE);
                    break;
            }
        }


        var mappings = FilterMappings.parseAndFilterMappingFile(mappingFile, allowedMappings, requireAll);
        FilterMappings.writeMappingFile(mappings, cmd.getOptionValue("output-file"));
    }
}
