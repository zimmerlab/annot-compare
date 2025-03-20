package com.github.zimmerlab.gtfcompare;

import com.github.kleinsamuel.gtfutils.GtfConfig;
import com.github.kleinsamuel.gtfutils.GtfConstants;
import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.GeneFeature;
import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import org.apache.commons.cli.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.boot.CommandLineRunner;
import org.springframework.context.annotation.Profile;
import org.springframework.stereotype.Service;

import java.io.File;
import java.util.List;
import java.util.Set;

@Profile("test")
@Service
public class TestRunner implements CommandLineRunner {

    private final static Logger LOG = LogManager.getLogger(TestRunner.class);

    public TestRunner() {

    }

    @Override
    public void run(String... args) throws Exception {

        Options o = new Options();
        o.addOption(Option.builder()
                .option("h")
                .longOpt("help")
                .desc("Print the help message")
                .build());
        o.addOption(Option.builder()
                .longOpt("gtf")
                .numberOfArgs(1)
                .required()
                .desc("Path to gtf file")
                .type(File.class)
                .build());

        CommandLineParser parser = new DefaultParser();

        CommandLine cmd = null;

        try {
            cmd = parser.parse(o, args);
        } catch (ParseException e) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Test", o, true);
            System.exit(1);
        }

        LOG.info("Running test");

        if (!cmd.hasOption("gtf")) {
            LOG.error("No gtf file specified");
            System.exit(1);
        }

        GtfFile gtfFile = new GtfFile(new File(cmd.getOptionValue("gtf")));

        // parse all contigs from the gtf file
//        gtfFile.parseAllContigs();

        // parse the next (first in this case) contig (chr) from the gtf file
        gtfFile.parseNextContig();

        LOG.info("parsed contig: {}", gtfFile.getParsedContig());

        Set<String> geneIds = gtfFile.getAllGeneFeatureIds();

        for (String geneId : geneIds) {

            LOG.info("Gene: {}", geneId);

            // every feature is a GtfFeature
            // GeneFeature is a subclass of GtfFeature that represents a gene
            GeneFeature geneFeature = gtfFile.getGeneFeature(geneId);

            LOG.info("Gene from feature: {}", geneFeature.getGeneId());

            // TranscriptFeature is a subclass of GtfFeature that represents a transcript
            List<TranscriptFeature> transcripts = geneFeature.getTranscripts();

            for (TranscriptFeature transcript : transcripts) {

                LOG.info("Transcript: {}", transcript.getTranscriptId());

                // GtfConfig contains the strings for each feature
                // each feature has some synonyms (e.g. exon, Exon, EXON, etc.) that are used for
                // parsing and then changed to the default which can be used for querying
                List<GtfFeature> exons = transcript.getFeatures(GtfConfig.TYPE_EXON_DEFAULT);

                for (GtfFeature exon : exons) {

                    // each attribute of the gtf can be queried by GtfConstants
                    // gtf allows multiple values per attribute key
                    String exonId = exon.getBaseData().getAttributes(GtfConstants.EXON_ID_ATTRIBUTE_KEY).get(0);

                    // the base attributes of the gtf feature (columns except attributes)
                    // are stored in the getBaseData() object
                    int exonStart = exon.getBaseData().getStart();
                    int exonEnd = exon.getBaseData().getEnd();

                    LOG.info("Exon: {} from {} to {}", exonId, exonStart, exonEnd);

                    break;
                }

                break;
            }

            break;
        }
    }
}
