package com.github.zimmerlab.gtfcompare.runner;

import com.github.kleinsamuel.gtfutils.GtfConfig;
import com.github.kleinsamuel.gtfutils.GtfFile;
import org.apache.commons.cli.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.boot.CommandLineRunner;
import org.springframework.context.annotation.Profile;
import org.springframework.stereotype.Service;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

@Profile("gtfStats")
@Service
public class GtfStatsRunner implements CommandLineRunner {
    private final static Logger logger = LogManager.getLogger(GtfStatsRunner.class);

    public void run(String... args) throws Exception {
        Options o = new Options();
        o.addOption(Option.builder().option("h").longOpt("help").desc("Print the help message").build());
        o.addOption(Option.builder().longOpt("inDir").hasArg().required().desc("Directory with one or more GTF files").type(File.class).build());
        o.addOption(Option.builder().longOpt("outDir").hasArg().required().desc("Directory to write per-GTF stats").type(File.class).build());
        o.addOption(Option.builder().longOpt("recursive").optionalArg(true).desc("Recursively scan subdirectories").build());

        CommandLineParser parser = new DefaultParser();
        CommandLine cmd;
        try {
            cmd = parser.parse(o, args);
        } catch (ParseException e) {
            new HelpFormatter().printHelp("gtfStats", o, true);
            System.exit(1);
            return;
        }

        var inDir = Path.of(cmd.getOptionValue("inDir"));
        var outDir = Path.of(cmd.getOptionValue("outDir"));
        var recursive = cmd.hasOption("recursive");

        if (!Files.isDirectory(inDir)) {
            logger.error("Input is not a directory: {}", inDir);
            System.exit(1);
        }
        Files.createDirectories(outDir);

        try (var stream = (recursive ? Files.walk(inDir) : Files.list(inDir))) {
            stream.filter(Files::isRegularFile).filter(p -> isGtfLike(p.getFileName().toString())).sorted().forEach(path -> {
                try {
                    String prefix = derivePrefix(path.getFileName().toString());
                    logger.info("Processing {} -> prefix {}", path, prefix);

                    var gtfFile = new GtfFile(path.toFile());
                    getStatsPerFile(gtfFile, outDir.toString(), prefix);

                } catch (Exception ex) {
                    logger.error("Failed processing {}: {}", path, ex.toString(), ex);
                }
            });
        }
    }


    static void getStatsPerFile(GtfFile gtfFile, String outputDir, String prefix) throws IOException {
        int numGenes = 0;
        int numTranscripts = 0;

        int numExons = 0;
        int numIntrons = 0;
        int numCDS = 0;
        int fiveUTR = 0;
        int threeUTR = 0;
        int numStartCodons = 0;
        int numStopCodons = 0;

        var geneBiotypeCount = new HashMap<String, Integer>();
        var transcriptBiotypeCount = new HashMap<String, Integer>();
        var transcriptsPerGene = new ArrayList<Integer>();
        var exonsPerTranscript = new ArrayList<Integer>();
        var intronsPerTranscript = new ArrayList<Integer>();

        var geneLengths = new ArrayList<Integer>();
        var transcriptLengths = new ArrayList<Integer>();
        var exonLengths = new ArrayList<Integer>();
        var intronLengths = new ArrayList<Integer>();
        var cdsLengths = new ArrayList<Integer>();
        var utr5Lengths = new ArrayList<Integer>();
        var utr3Lengths = new ArrayList<Integer>();

        try {
            while (true) {
                gtfFile.parseNextContig();

                for (var geneId : gtfFile.getAllGeneFeatureIds()) {
                    var gene = gtfFile.getGeneFeature(geneId);
                    numGenes++;
                    transcriptsPerGene.add(gene.getTranscripts().size());

                    var geneBaseData = gene.getBaseData();
                    geneLengths.add(geneBaseData.getEnd() - geneBaseData.getStart() + 1);

                    var geneBiotype = geneBaseData.getAttributes("gene_biotype") == null ? geneBaseData.getSource() : geneBaseData.getAttributes("gene_biotype").getFirst();
                    geneBiotypeCount.merge(geneBiotype, 1, Integer::sum);

                    for (var transcript : gene.getTranscripts()) {
                        numTranscripts++;
                        var transcriptBaseData = transcript.getBaseData();
                        var tStart = transcriptBaseData.getStart();
                        var tEnd = transcriptBaseData.getEnd();

                        if(tStart == -1 || tEnd == -1) {
                            tStart = transcript.getFeatures().stream()
                                    .mapToInt(query -> query.getBaseData().getStart())
                                    .min()
                                    .orElse(Integer.MAX_VALUE);

                            tEnd = transcript.getFeatures().stream()
                                    .mapToInt(query -> query.getBaseData().getEnd())
                                    .max()
                                    .orElse(Integer.MIN_VALUE);
                        }
                        transcriptLengths.add(tEnd - tStart + 1);
                        var transcriptBiotype = transcriptBaseData.getAttributes("transcript_biotype") == null ? transcriptBaseData.getSource() : transcriptBaseData.getAttributes("transcript_biotype").getFirst();
                        transcriptBiotypeCount.merge(transcriptBiotype, 1, Integer::sum);

                        var featureCounts = new HashMap<String, Integer>();
                        for (var transcriptFeature : transcript.getFeatures()) {
                            var transcriptFeatureBaseData = transcriptFeature.getBaseData();
                            int len = transcriptFeatureBaseData.getEnd() - transcriptFeatureBaseData.getStart() + 1;
                            var type = GtfConfig.getDefault(transcriptFeatureBaseData.getType());
                            featureCounts.merge(type, 1, Integer::sum);

                            switch (type) {
                                case "exon":
                                    numExons++;
                                    exonLengths.add(len);
                                    break;
                                case "intron":
                                    numIntrons++;
                                    intronLengths.add(len);
                                    break;
                                case "CDS":
                                    numCDS++;
                                    cdsLengths.add(len);
                                    break;
                                case "five_prime_utr":
                                    fiveUTR++;
                                    utr5Lengths.add(len);
                                    break;
                                case "three_prime_utr":
                                    threeUTR++;
                                    utr3Lengths.add(len);
                                    break;
                                case "start_codon":
                                    numStartCodons++;
                                    break;
                                case "stop_codon":
                                    numStopCodons++;
                                    break;
                                default:
                                    break;
                            }
                        }

                        exonsPerTranscript.add(featureCounts.getOrDefault("exon", 0));
                        intronsPerTranscript.add(featureCounts.getOrDefault("intron", 0));
                    }
                }
            }
        } catch (java.text.ParseException e) {
            logger.info("Calculations finished");
        } catch (Exception e) {
            logger.error("Program failed:", e);
        }

        var outDir = Path.of(outputDir);
        Files.createDirectories(outDir);

        var globalPath = outDir.resolve("global_counts").resolve(prefix + "_global_counts.tsv");
        Files.createDirectories(globalPath.getParent());
        try (var pw = new PrintWriter(Files.newBufferedWriter(globalPath))) {
            pw.println("metric\tvalue");
            pw.printf("genes\t%d%n", numGenes);
            pw.printf("transcripts\t%d%n", numTranscripts);
            pw.printf("exons\t%d%n", numExons);
            pw.printf("introns\t%d%n", numIntrons);
            pw.printf("cds\t%d%n", numCDS);
            pw.printf("five_prime_utr\t%d%n", fiveUTR);
            pw.printf("three_prime_utr\t%d%n", threeUTR);
            pw.printf("start_codons\t%d%n", numStartCodons);
            pw.printf("stop_codons\t%d%n", numStopCodons);
        }

        writeMapTsv(outDir.resolve("gene_biotypes").resolve(prefix + "_gene_biotypes.tsv"), "biotype", "count", geneBiotypeCount);
        writeMapTsv(outDir.resolve("transcript_biotypes").resolve(prefix + "_transcript_biotypes.tsv"), "biotype", "count", transcriptBiotypeCount);

        writeListTsv(outDir.resolve("transcripts_per_gene").resolve(prefix + "_transcripts_per_gene.tsv"), "transcripts_per_gene", transcriptsPerGene);
        writeListTsv(outDir.resolve("exons_per_transcript").resolve(prefix + "_exons_per_transcript.tsv"), "exons_per_transcript", exonsPerTranscript);
        writeListTsv(outDir.resolve("introns_per_transcript").resolve(prefix + "_introns_per_transcript.tsv"), "introns_per_transcript", intronsPerTranscript);

        writeListTsv(outDir.resolve("exon_lengths").resolve(prefix + "_exon_lengths.tsv"), "exon_length", exonLengths);
        writeListTsv(outDir.resolve("intron_lengths").resolve(prefix + "_intron_lengths.tsv"), "intron_length", intronLengths);
        writeListTsv(outDir.resolve("cds_lengths").resolve(prefix + "_cds_lengths.tsv"), "cds_length", cdsLengths);
        writeListTsv(outDir.resolve("utr5_lengths").resolve(prefix + "_utr5_lengths.tsv"), "utr5_length", utr5Lengths);
        writeListTsv(outDir.resolve("utr3_lengths").resolve(prefix + "_utr3_lengths.tsv"), "utr3_length", utr3Lengths);
        writeListTsv(outDir.resolve("gene_lengths").resolve(prefix + "_gene_lengths.tsv"), "gene_length", geneLengths);
        writeListTsv(outDir.resolve("transcript_lengths").resolve(prefix + "_transcript_lengths.tsv"), "transcript_length", transcriptLengths);

    }

    private static void writeMapTsv(Path path, String keyHeader, String valHeader, Map<String, Integer> map) throws IOException {
        Files.createDirectories(path.getParent());
        try (var pw = new PrintWriter(Files.newBufferedWriter(path))) {
            pw.println(keyHeader + "\t" + valHeader);
            map.entrySet().stream().sorted(Map.Entry.comparingByKey()).forEach(e -> pw.printf("%s\t%d%n", e.getKey(), e.getValue()));
        }
    }

    private static void writeListTsv(Path path, String header, List<Integer> values) throws IOException {
        Files.createDirectories(path.getParent());
        try (var pw = new PrintWriter(Files.newBufferedWriter(path))) {
            pw.println(header);
            for (int v : values) pw.println(v);
        }
    }

    private static boolean isGtfLike(String name) {
        String n = name.toLowerCase();
        return n.endsWith(".gtf");
    }

    private static String derivePrefix(String filename) {
        String n = filename;
        var ext = ".gtf";
        if (n.toLowerCase().endsWith(ext)) {
            n = n.substring(0, n.length() - ext.length());
        }
        return n;
    }

}
