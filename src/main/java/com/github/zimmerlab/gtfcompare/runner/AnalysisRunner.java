package com.github.zimmerlab.gtfcompare.runner;

import com.fasterxml.jackson.databind.MapperFeature;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.AnnotComparator;
import com.github.zimmerlab.gtfcompare.compare.ComparisonConfig;
import com.github.zimmerlab.gtfcompare.compare.ComparisonConfigBuilder;
import com.github.zimmerlab.gtfcompare.mapping.Minimap2Bundler;
import com.github.zimmerlab.gtfcompare.mapping.Minimap2Validator;
import com.github.zimmerlab.gtfcompare.mapping.OverlappingTranscripts;
import com.github.zimmerlab.gtfcompare.model.MappingResult;
import com.github.zimmerlab.gtfcompare.model.TranscriptPair;
import com.github.zimmerlab.gtfcompare.model.config.ConfigJSON;
import com.github.zimmerlab.gtfcompare.model.config.FeatureConfig;
import com.github.zimmerlab.gtfcompare.parser.FidxParser;
import com.github.zimmerlab.gtfcompare.utils.*;
import org.apache.commons.cli.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.boot.CommandLineRunner;
import org.springframework.context.annotation.Profile;
import org.springframework.stereotype.Service;
import org.springframework.util.StopWatch;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;

@Profile("analysis")
@Service
public class AnalysisRunner implements CommandLineRunner {
    private final static Logger logger = LogManager.getLogger(AnalysisRunner.class);
    private final static List<StopWatch> stopWatches = Constants.STOP_WATCHES;

    @Override
    public void run(String... args) throws Exception {
        Options o = new Options();
        o.addOption(Option.builder().option("h").longOpt("help").desc("Print the help message").build());
        o.addOption(Option.builder().longOpt("gtf").numberOfArgs(1).required().desc("Path to gtf file").type(File.class).build());

        o.addOption(Option.builder().longOpt("fasta").numberOfArgs(1).required().desc("Path to fasta file").type(File.class).build());

        o.addOption(Option.builder().longOpt("fidx").numberOfArgs(1).required().desc("Path to fasta index file").type(File.class).build());

        o.addOption(Option.builder().longOpt("gtf2").numberOfArgs(1).required().desc("Path to gtf file").type(File.class).build());

        o.addOption(Option.builder().longOpt("fasta2").numberOfArgs(1).required().desc("Path to fasta file").type(File.class).build());

        o.addOption(Option.builder().longOpt("fidx2").numberOfArgs(1).required().desc("Path to fasta index file").type(File.class).build());

        o.addOption(Option.builder().longOpt("o").numberOfArgs(1).required().desc("Path to output file").type(File.class).build());

        o.addOption(Option.builder().longOpt("config").numberOfArgs(1).required().desc("Path to config file").type(File.class).build());

        /*o.addOption(Option.builder()
                .longOpt("gene-mapping")
                .numberOfArgs(1)
                .required()
                .desc("Path to mapping file")
                .type(File.class)
                .build());*/

        CommandLineParser parser = new DefaultParser();

        CommandLine cmd = null;

        try {
            cmd = parser.parse(o, args);
        } catch (ParseException e) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Test", o, true);
            System.exit(1);
        }

        logger.info("Running test");

        if (!cmd.hasOption("gtf")) {
            logger.error("No gtf file specified");
            System.exit(1);
        }

        if (!cmd.hasOption("fasta")) {
            logger.error("No fasta file specified");
            System.exit(1);
        }

        if (!cmd.hasOption("fidx")) {
            logger.error("No fidx file specified");
            System.exit(1);
        }

        if (!cmd.hasOption("config")) {
            logger.error("No config file specified");
            System.exit(1);
        }

        if (!cmd.hasOption("o")) {
            logger.error("No output path specified");
            System.exit(1);
        }

        var useLiftoff = false;

        if (useLiftoff) {
            String userHome = System.getProperty("user.home");
            String liftoffExec = Paths.get(userHome, "miniconda3", "bin", "liftoff").toString();

            ProcessBuilder pb = new ProcessBuilder(liftoffExec, "-g", cmd.getOptionValue("gtf"), cmd.getOptionValue("fasta"), cmd.getOptionValue("fasta2"), "-o", "output/lifted.gtf");

            pb.redirectErrorStream(true);

            try {
                Process process = pb.start();

                try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
                    String line;
                    while ((line = reader.readLine()) != null) {
                        System.out.println(line);
                    }
                }

                int exitCode = process.waitFor();
                if (exitCode == 0) {
                    System.out.println("Liftoff successfully executed.");
                } else {
                    System.err.println("Liftoff failed with exit code " + exitCode + ".");
                }

            } catch (IOException | InterruptedException e) {
                e.printStackTrace();
            }
        }

        var fidxEntries = FidxParser.parse(cmd.getOptionValue("fidx"));
        var fidx2Entries = FidxParser.parse(cmd.getOptionValue("fidx2"));

        GtfFile gtfFile = new GtfFile(new File(cmd.getOptionValue("gtf")));
        GtfFile gtfFile2 = new GtfFile(new File(cmd.getOptionValue("gtf2")));

        var targetSequenceExtractor = new GenomeSequenceExtractor(new File(cmd.getOptionValue("fasta")), fidxEntries);
        var querySequenceExtractor = new GenomeSequenceExtractor(new File(cmd.getOptionValue("fasta2")), fidx2Entries);

        var queryUnmappedFilePath = Path.of("output/unmapped_queries.txt");
        var targetUnmappedFilePath = Path.of("output/unmapped_targets.txt");
        try (BufferedWriter writer = Files.newBufferedWriter(queryUnmappedFilePath,
                StandardOpenOption.CREATE,
                StandardOpenOption.TRUNCATE_EXISTING)) {
            writer.write("");
        }
        try (BufferedWriter writer = Files.newBufferedWriter(targetUnmappedFilePath,
                StandardOpenOption.CREATE,
                StandardOpenOption.TRUNCATE_EXISTING)) {
            writer.write("");
        }

        final String HEADER = String.join("\t", "impact", "contig", "targetGeneId", "queryGeneId", "targetBioType", "queryBiotype", "featureType", "difference", "targetTranscriptId", "queryTranscriptId", "targetTranscriptBiotype", "queryTranscriptBiotype", "targetFeatureStart", "queryFeatureStart", "targetFeatureStop", "queryFeatureStop", "targetStrand", "queryStrand");

        try (BufferedWriter writer = Files.newBufferedWriter(Path.of(cmd.getOptionValue("o")),
                StandardOpenOption.CREATE,
                StandardOpenOption.TRUNCATE_EXISTING)) {
            writer.write(HEADER);
            writer.newLine();
        }

        try (BufferedWriter writer = Files.newBufferedWriter(Path.of(cmd.getOptionValue("o") + ".minimap2"),
                StandardOpenOption.CREATE,
                StandardOpenOption.TRUNCATE_EXISTING)) {
            writer.write(HEADER);
            writer.newLine();
        }


        try {
            while (true) {
                gtfFile.parseNextContig();
                gtfFile2.parseNextContig();
                if (!gtfFile.getParsedContig().equals(gtfFile2.getParsedContig())) {
                    throw new Exception("Contigs do not match");
                }
                runProgrammePerContig(gtfFile, gtfFile2, targetSequenceExtractor, querySequenceExtractor, cmd, targetUnmappedFilePath, queryUnmappedFilePath);
            }
        } catch (java.text.ParseException e) {
            logger.info("Program finished");
        } catch (Exception e) {
            logger.error("Program failed", e);
        }

    }

    private static void runProgrammePerContig(GtfFile gtfFile, GtfFile gtfFile2, GenomeSequenceExtractor targetSequenceExtractor, GenomeSequenceExtractor querySequenceExtractor, CommandLine cmd, Path targetUnmappedPath, Path queryUnmappedPath) throws IOException, InterruptedException {
        var currentContig = gtfFile.getParsedContig();
        logger.info("Current contig: " + currentContig);
        var loci = OverlappingTranscripts.map(gtfFile, gtfFile2);



        var workDir = Files.createTempDirectory("mm2-");
        workDir.toFile().deleteOnExit();

        var minimapPath = Minimap2Bundler.extractMinimap2();
        logger.debug("Starting Minimap");
        var minimap2Result = Minimap2Validator.validateWithMinimap2(loci.getUnmappedQueries(), loci.getUnmappedTargets(), targetSequenceExtractor, querySequenceExtractor, workDir, minimapPath, 8);
        try (BufferedWriter writer = Files.newBufferedWriter(queryUnmappedPath, StandardOpenOption.CREATE, StandardOpenOption.APPEND)) {
            for (TranscriptFeature tf : minimap2Result.getUnmappedQueries()) {
                writer.write(tf.getBaseData().getContig() + "\t" + tf.getTranscriptId());
                writer.newLine();
            }
        }

        try (BufferedWriter writer = Files.newBufferedWriter(targetUnmappedPath, StandardOpenOption.CREATE, StandardOpenOption.APPEND)) {
            for (TranscriptFeature tf : minimap2Result.getUnmappedTargets()) {
                writer.write(tf.getBaseData().getContig() + "\t" + tf.getTranscriptId());
                writer.newLine();
            }
        }

        var configPath = cmd.getOptionValue("config");

        var objectMapper = new ObjectMapper();
        objectMapper.configure(MapperFeature.ACCEPT_CASE_INSENSITIVE_ENUMS, true);
        ConfigJSON jsonConfig = null;
        try {
            jsonConfig = objectMapper.readValue(new File(configPath), ConfigJSON.class);
        } catch (Exception e) {
            logger.error("Error reading config file: {}", e.getMessage());
            System.exit(1);
        }

        var config = getComparisonConfig(jsonConfig);

        logger.debug("Starting with normal comparison");
        var annotCompareSafe = new AnnotComparator(gtfFile, gtfFile2, targetSequenceExtractor, querySequenceExtractor, config, cmd.getOptionValue("o"), loci.getMapping());
        annotCompareSafe.compare();

        logger.debug("Starting with minimap comparison");
        var annotCompareMm2 = new AnnotComparator(gtfFile, gtfFile2, targetSequenceExtractor, querySequenceExtractor, config, cmd.getOptionValue("o") + ".minimap2", minimap2Result.getMapping());
        annotCompareMm2.compare();
    }

    private static ComparisonConfig getComparisonConfig(ConfigJSON jsonConfig) {
        var configBuilder = new ComparisonConfigBuilder();
        configBuilder.setAllowedGeneBiotypes(jsonConfig.getGeneBiotypes().getOrDefault("allowed", new ArrayList<>()));
        var configFeatures = jsonConfig.getFeatures();

        // FEATURE COMPARATORS
        enableFeatureWithThreshold(configBuilder, Constants.LENGTH_COMPARATOR_NAME, configFeatures);
        enableFeatureWithThreshold(configBuilder, Constants.START_COMPARATOR_NAME, configFeatures);
        enableFeatureWithThreshold(configBuilder, Constants.STOP_COMPARATOR_NAME, configFeatures);

        enableFeature(configBuilder, Constants.SEQUENCE_COMPARATOR_NAME, configFeatures);
        enableFeature(configBuilder, Constants.SAME_PROTEIN_COMPARATOR_NAME, configFeatures);


        // TRANSCRIPT COMPARATORS

        enableFeatureWithThreshold(configBuilder, Constants.TRANSCRIPT_LENGTH_COMPARATOR_NAME, configFeatures);
        enableFeatureWithThreshold(configBuilder, Constants.TRANSCRIPT_START_COMPARATOR_NAME, configFeatures);
        enableFeatureWithThreshold(configBuilder, Constants.TRANSCRIPT_STOP_COMPARATOR_NAME, configFeatures);
        enableFeatureWithThreshold(configBuilder, Constants.TRANSCRIPT_STRAND_COMPARATOR_NAME, configFeatures);
        enableFeatureWithThreshold(configBuilder, Constants.TRANSCRIPT_BIOTYPE_COMPARATOR_NAME, configFeatures);


        // GENE COMPARATORS

        enableFeatureWithThreshold(configBuilder, Constants.GENE_LENGTH_COMPARATOR_NAME, configFeatures);
        enableFeatureWithThreshold(configBuilder, Constants.START_COMPARATOR_NAME, configFeatures);
        enableFeatureWithThreshold(configBuilder, Constants.STOP_COMPARATOR_NAME, configFeatures);

        enableFeature(configBuilder, Constants.GENE_STRAND_COMPARATOR_NAME, configFeatures);
        enableFeature(configBuilder, Constants.GENE_CONTIG_COMPARATOR_NAME, configFeatures);

        var transcriptFeatures = jsonConfig.getTranscriptFeatures();

        for (var transcriptFeature : Constants.FEATURE_TYPES) {
            var feature = transcriptFeatures.get(transcriptFeature);
            if (feature != null && feature.isEnabled()) {
                configBuilder.enableTranscriptFeatures(transcriptFeature);
                var impactLvl = feature.getImpactLevel();
                if(impactLvl != null){
                    configBuilder.setImpactLevels(transcriptFeature, impactLvl);

                }
                /*var th = feature.getThreshold();
                if (th != null) {
                    configBuilder.setThreshold(transcriptFeature, th);
                }*/
            }
        }

        return configBuilder.build();
    }

    private static void enableFeature(ComparisonConfigBuilder configBuilder, String featureName, Map<String, FeatureConfig> featureConfig) {
        var feature = featureConfig.get(featureName);
        if (feature != null && feature.isEnabled()) {
            configBuilder.enableFeature(featureName);
        }

        if(feature.getImpactLevel() == null) return;
        configBuilder.setImpactLevels(featureName, feature.getImpactLevel() == null ? null : feature.getImpactLevel());
    }

    private static void enableFeatureWithThreshold(ComparisonConfigBuilder configBuilder, String featureName, Map<String, FeatureConfig> featureConfig) {
        var feature = featureConfig.get(featureName);
        if (feature != null && feature.isEnabled()) {
            configBuilder.enableFeature(featureName);
            var th = feature.getThreshold();
            if (th != null) {
                configBuilder.setThreshold(featureName, th);
            }
        }

        if(feature.getImpactLevel() == null) return;
        configBuilder.setImpactLevels(featureName, feature.getImpactLevel());
    }
    private static void WriteUnmappedAfterOverlaps(MappingResult<TranscriptPair, TranscriptFeature> loci) {
        try (BufferedWriter writer = Files.newBufferedWriter(Path.of("unmapped.tsv"))) {
            writer.write("source\ttranscript_id\tfeature\tstart\tstop\n");

            for (var unmapped : loci.getUnmappedQueries()) {
                var baseData = unmapped.getBaseData();

                writer.write("query\t");
                var transcriptId = baseData.getAttributes("transcript_id").get(0);
                writer.write(transcriptId + "\t" + baseData.getType() + "\t" + baseData.getStart() + "\t" + baseData.getEnd() + "\n");

                for (var feature : unmapped.getFeatures()) {
                    var featureBaseData = feature.getBaseData();
                    writer.write("query\t");
                    writer.write(transcriptId + "\t" + featureBaseData.getType() + "\t" + featureBaseData.getStart() + "\t" + featureBaseData.getEnd() + "\n");
                }
            }

            for (var unmapped : loci.getUnmappedTargets()) {
                var baseData = unmapped.getBaseData();

                writer.write("target\t");
                var transcriptId = baseData.getAttributes("transcript_id").get(0);
                writer.write(transcriptId + "\t" + baseData.getType() + "\t" + baseData.getStart() + "\t" + baseData.getEnd() + "\n");

                for (var feature : unmapped.getFeatures()) {
                    var featureBaseData = feature.getBaseData();
                    writer.write("target\t");
                    writer.write(transcriptId + "\t" + featureBaseData.getType() + "\t" + featureBaseData.getStart() + "\t" + featureBaseData.getEnd() + "\n");
                }
            }

        } catch (IOException e) {
            // TODO logging
            throw new RuntimeException(e);
        }
    }

}
