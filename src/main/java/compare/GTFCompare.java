package compare;

import com.github.kleinsamuel.gtfutils.feature.GeneFeature;
import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.GenomeSequenceExtractor;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.*;

public class GTFCompare {
    public static final Logger LOG = LoggerFactory.getLogger(GTFCompare.class);

    public static void compare(GeneFeature gene1, GeneFeature gene2, GenomeSequenceExtractor sequenceExtractor1, GenomeSequenceExtractor sequenceExtractor2) {
        var baseData1 = gene1.getBaseData();
        var baseData2 = gene2.getBaseData();

        if (baseData1.getStart() != baseData2.getStart()) {
            LOG.info("Start position not the same");
        }

        if (baseData2.getEnd() != baseData2.getEnd()) {
            LOG.info("Stop position not the same");
        }

        // TODO für alle oder nur gen?
        if(baseData1.isForwardStrand() != baseData2.isForwardStrand()){
            LOG.info("Strand changed from {} to {}", baseData1.isForwardStrand(), baseData2.isForwardStrand());
        }

        try {
            var seq1 = sequenceExtractor1.getSequence(baseData1.getContig(), baseData1.getStart(), baseData1.getEnd());
            var seq2 = sequenceExtractor2.getSequence(baseData2.getContig(), baseData2.getStart(), baseData2.getEnd());

            if(!seq1.equals(seq2)){
               LOG.info("Gene sequence not the same");
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }


        compareTranscripts(gene1, gene2);

    }

    private static void compareTranscripts(GeneFeature g1, GeneFeature g2) {

        var transcriptsMap1 = new HashMap<String, TranscriptFeature>();
        var transcriptsMap2 = new HashMap<String, TranscriptFeature>();

        for (var transcript : g1.getTranscripts()) {
            transcriptsMap1.computeIfAbsent(transcript.getTranscriptId(), key -> transcript);
        }

        for (var transcript : g2.getTranscripts()) {
            transcriptsMap2.computeIfAbsent(transcript.getTranscriptId(), key -> transcript);
        }

        for (String transcriptId : transcriptsMap1.keySet()) {
            if (!transcriptsMap2.containsKey(transcriptId)) {
                LOG.info("Transcript {} missing in gene2.", transcriptId);
            } else {
                TranscriptFeature t1 = transcriptsMap1.get(transcriptId);
                TranscriptFeature t2 = transcriptsMap2.get(transcriptId);

                var baseData1 = t1.getBaseData();
                var baseData2 = t2.getBaseData();

                if (baseData1.getStart() != baseData2.getStart()) {
                    LOG.info("Start position not the same");
                }

                if (baseData2.getEnd() != baseData2.getEnd()) {
                    LOG.info("Stop position not the same");
                }

                compareFeatures(t1, t2);
            }
        }

        for (String transcriptId : transcriptsMap2.keySet()) {
            if (!transcriptsMap1.containsKey(transcriptId)) {
                LOG.info("Transcript {} missing in gene1.", transcriptId);
            }
        }
    }

    private static void compareFeatures(TranscriptFeature t1, TranscriptFeature t2) {
        var featureMap1 = new HashMap<String, List<GtfFeature>>();
        var featureMap2 = new HashMap<String, List<GtfFeature>>();

        for (var feature : t1.getFeatures()) {
            featureMap1.computeIfAbsent(feature.getBaseData().getType(), key -> new ArrayList<>()).add(feature);
        }

        for (var feature : t2.getFeatures()) {
            featureMap2.computeIfAbsent(feature.getBaseData().getType(), key -> new ArrayList<>()).add(feature);
        }

        compareFeaturePositions(t1, t2, featureMap1, featureMap2);

    }

    private static void compareFeaturePositions(TranscriptFeature t1, TranscriptFeature t2, Map<String, List<GtfFeature>> featureMap1, Map<String, List<GtfFeature>> featureMap2) {
        for (var featureKey : featureMap1.keySet()) {
            if (!featureMap2.containsKey(featureKey)) {
                LOG.info("{} not in {}", featureKey, t2.getTranscriptId());
                continue;
            }
            List<GtfFeature> sortedExons1 = new ArrayList<>(featureMap1.get(featureKey));
            List<GtfFeature> sortedExons2 = new ArrayList<>(featureMap2.get(featureKey));

            sortedExons1.sort(
                    Comparator.comparingInt((GtfFeature f) -> f.getBaseData().getStart())
                            .thenComparingInt(f -> f.getBaseData().getEnd())
            );
            sortedExons2.sort(
                    Comparator.comparingInt((GtfFeature f) -> f.getBaseData().getStart())
                            .thenComparingInt(f -> f.getBaseData().getEnd())
            );

            int i = 0, j = 0;
            while (i < sortedExons1.size() && j < sortedExons2.size()) {
                GtfFeature exon1 = sortedExons1.get(i);
                GtfFeature exon2 = sortedExons2.get(j);

                int start1 = exon1.getBaseData().getStart();
                int start2 = exon2.getBaseData().getStart();
                int end1 = exon1.getBaseData().getEnd();
                int end2 = exon2.getBaseData().getEnd();

                // Case: Both start and end are identical.
                if (start1 == start2 && end1 == end2) {
                    i++;
                    j++;

                }
                // Case: Different end position
                else if (start1 == start2) {
                    LOG.info("{} end changed: list 1: {}-{}, list 2: {}-{}", featureKey, start1, end1, start2, end2);
                    i++;
                    j++;
                }
                // Case: Different start position
                else if (end1 == end2) {
                    LOG.info("{} start changed: list 1: {}-{}, list 2: {}-{}", featureKey, start1, end1, start2, end2);
                    i++;
                    j++;
                }
                // Case: exon1 has an earlier start, indicating it is missing in list2.
                else if (start1 < start2) {
                    LOG.info("{} in {} missing in list 2: {}-{}", featureKey, t1.getTranscriptId(), start1, end1);
                    i++;
                }
                // Otherwise: exon2 must have an earlier start, so report it as missing in list1.
                else {
                    LOG.info("{} in {} missing in list 1: {}-{}", featureKey, t2.getTranscriptId(), start2, end2);
                    j++;
                }
            }

            // Log extra features if any remain.
            while (i < sortedExons1.size()) {
                GtfFeature exon1 = sortedExons1.get(i);
                LOG.info("Extra {} in {} in list 1: {}-{}", featureKey, t1.getTranscriptId(), exon1.getBaseData().getStart(), exon1.getBaseData().getEnd());
                i++;
            }
            while (j < sortedExons2.size()) {
                GtfFeature exon2 = sortedExons2.get(j);
                LOG.info("Extra {} in {} in list 2: {}-{}", featureKey, t2.getTranscriptId(), exon2.getBaseData().getStart(), exon2.getBaseData().getEnd());
                j++;
            }
        }
    }

}
