package com.github.zimmerlab.gtfcompare.runner;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.GtfBaseData;
import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
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
import java.nio.file.StandardOpenOption;
import java.util.List;
import java.util.stream.Collectors;

@Profile("addMetaFeatures")
@Service
public class AddMetaFeaturesToGTFRunner implements CommandLineRunner {

    private final static Logger logger = LogManager.getLogger(AddMetaFeaturesToGTFRunner.class);

    @Override
    public void run(String... args) throws Exception {

        Options o = new Options();
        o.addOption(Option.builder().option("h").longOpt("help").desc("Print the help message").build());
        o.addOption(Option.builder().longOpt("gtf").numberOfArgs(1).required().desc("Path to gtf file").type(File.class).build());
        o.addOption(Option.builder().longOpt("o").numberOfArgs(1).required().desc("Path to gtf file").type(File.class).build());
        CommandLineParser parser = new DefaultParser();

        CommandLine cmd = null;

        try {
            cmd = parser.parse(o, args);
        } catch (ParseException e) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Test", o, true);
            System.exit(1);
        }

        if (!cmd.hasOption("gtf")) {
            logger.error("No gtf file specified");
            System.exit(1);
        }

        if (!cmd.hasOption("o")) {
            logger.error("No output file specified");
            System.exit(1);
        }

        var gtfFile = new GtfFile(new File(cmd.getOptionValue("gtf")));

        try (BufferedWriter writer = Files.newBufferedWriter(Path.of(cmd.getOptionValue("o")))) {

            while (true) {
                gtfFile.parseNextContig();

                var geneIds = gtfFile.getAllGeneFeatureIds();

                for (var geneId : geneIds) {
                    var gene = gtfFile.getGeneFeature(geneId);
                    var geneBaseData = gene.getBaseData();

                    var line = getGtfLine(geneBaseData, null);
                    writer.write(String.join("\t", line) + "\n");

                    for (var transcripts : gene.getTranscripts()) {
                        var transcriptBaseData = transcripts.getBaseData();
                        var line2 = getGtfLine(transcriptBaseData, transcripts.getFeatures());
                        writer.write(String.join("\t", line2) + "\n");

                        for(var transcriptFeature: transcripts.getFeatures()){
                            var transcriptFeatureBaseData = transcriptFeature.getBaseData();
                            if(transcriptFeatureBaseData.getType().equals("intron")){
                                continue;
                            }
                            var featureLine = getGtfLine(transcriptFeatureBaseData, null);
                            writer.write(String.join("\t", featureLine) + "\n");
                        }
                    }
                }
            }
        } catch (Exception e) {

        }
    }

    public static String[] getGtfLine(GtfBaseData baseData, List<GtfFeature> subfeatures) {

        var start = baseData.getStart();
        var end = baseData.getEnd();
        var type =  baseData.getType();

        if(type.equals("transcript") && (start == -1 || end == -1)){
            start = subfeatures.stream().mapToInt(target -> target.getBaseData().getStart())
                    .min()
                    .orElse(Integer.MAX_VALUE);

            end = subfeatures.stream().mapToInt(target -> target.getBaseData().getEnd())
                    .max()
                    .orElse(Integer.MIN_VALUE);
        }

        return new String[]{
                baseData.getContig(),
                baseData.getSource(),
                type,
                String.valueOf(start),
                String.valueOf(end),
                baseData.getScore() == null ? "." : String.valueOf(baseData.getScore()),
                baseData.isForwardStrand() ? "+" : "-",
                baseData.getFrame() == null ? "." : String.valueOf(baseData.getFrame()),
                baseData.getAttributes().entrySet().stream().map(entry -> entry.getKey() + " \"" + (!entry.getValue().isEmpty() ? entry.getValue().getFirst() + "\"" : "")).collect(Collectors.joining(";")),
        };
    }
}
