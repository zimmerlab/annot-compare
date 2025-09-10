package com.github.zimmerlab.gtfcompare;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.boot.WebApplicationType;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.boot.builder.SpringApplicationBuilder;

import java.util.Properties;
import java.util.TreeMap;

@SpringBootApplication
public class Application {

    private final static Logger LOG  = LogManager.getLogger(Application.class);

    private static void printModeUsage(TreeMap<String, String> modeDescriptions) {
        System.out.println("Usage: <mode> [options]");
        System.out.println("Modes:");
        for (String mode : modeDescriptions.keySet()) {
            System.out.printf("\t%s - %s\n", mode, modeDescriptions.get(mode));
        }
    }

    private static TreeMap<String, String> getModeDescriptions() {
        TreeMap<String, String> modeDescriptions = new TreeMap<>();

        modeDescriptions.put("test", "Test mode");
        modeDescriptions.put("analysis", "Analysis mode");
        modeDescriptions.put("firstAnalysis", "First Analysis mode");
        modeDescriptions.put("seqExtractor", "Sequence extractor mode");
        modeDescriptions.put("addMetaFeatures", "Add Meta Features mode");
        modeDescriptions.put("benchmark", "Benchmark mode");

        return modeDescriptions;
    }

    public static void main(String[] args) {

        TreeMap<String, String> modeDescriptions = getModeDescriptions();
        if (args.length == 0) {
            LOG.error("No mode specified\n");
            printModeUsage(modeDescriptions);
            System.exit(1);
        }

        SpringApplicationBuilder builder = new SpringApplicationBuilder(Application.class);
        Properties properties = new Properties();

        String mode = args[0];

        switch (mode) {
            case "test":
                builder.profiles("test");
                builder.web(WebApplicationType.NONE);
                break;
            case "analysis":
                builder.profiles("analysis");
                builder.web(WebApplicationType.NONE);
                break;
            case "firstAnalysis":
                builder.profiles("firstAnalysis");
                builder.web(WebApplicationType.NONE);
                break;
            case "seqExtractor":
                builder.profiles("seqExtractor");
                builder.web(WebApplicationType.NONE);
                break;
            case "addMetaFeatures":
                builder.profiles("addMetaFeatures");
                builder.web(WebApplicationType.NONE);
                break;
            case "benchmark":
                builder.profiles("benchmark");
                builder.web(WebApplicationType.NONE);
                break;
            default:
                LOG.error("Unknown mode: {}", mode);
                printModeUsage(modeDescriptions);
                System.exit(1);
        }

        builder.properties(properties);
        builder.run(args);
    }

}
