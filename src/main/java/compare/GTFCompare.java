package compare;

import com.github.kleinsamuel.gtfutils.feature.GeneFeature;
import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;

public class GTFCompare {
    public static final Logger LOG = LoggerFactory.getLogger(GTFCompare.class);
    public static void comparePosition(GeneFeature gene1, GeneFeature gene2){
        var baseData1 = gene1.getBaseData();
        var baseData2 = gene2.getBaseData();

        if(baseData1.getStart() != baseData2.getStart()){
            LOG.info("Start position not the same");
        }

        if(baseData2.getEnd() != baseData2.getEnd()){
            LOG.info("Stop position not the same");
        }

        compareTranscripts(gene1, gene2);

    }

    private static void compareTranscripts(GeneFeature g1, GeneFeature g2){

        var transcriptsMap1 = new HashMap<String, TranscriptFeature>();
        var transcriptsMap2 = new HashMap<String, TranscriptFeature>();

        for(var transcript : g1.getTranscripts()){
            transcriptsMap1.computeIfAbsent(transcript.getTranscriptId(), key -> transcript);
        }

        for(var transcript : g2.getTranscripts()){
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

                if(baseData1.getStart() != baseData2.getStart()){
                    LOG.info("Start position not the same");
                }

                if(baseData2.getEnd() != baseData2.getEnd()){
                    LOG.info("Stop position not the same");
                }

                compareFeatures(t1, t2);
            }
        }
    }

    private static void compareFeatures(TranscriptFeature t1, TranscriptFeature t2){
        var featureMap1 = new HashMap<String, List<GtfFeature>>();
        var featureMap2 = new HashMap<String, List<GtfFeature>>();

        for(var feature : t1.getFeatures()){
            featureMap1.computeIfAbsent(feature.getBaseData().getType(), key -> new ArrayList<>()).add(feature);
        }

        for(var feature : t2.getFeatures()){
            featureMap2.computeIfAbsent(feature.getBaseData().getType(), key -> new ArrayList<>()).add(feature);
        }

        for(var featureKey : featureMap1.keySet()){
            if(!featureMap2.containsKey(featureKey)){
                LOG.info(featureKey + " not in featuremap2");
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

                if (start1 == start2 && end1 == end2) {
                    i++;
                    j++;
                }
                // Fall: exon1 hat einen früheren Start, also fehlt dieses Exon vermutlich in list2
                else if (start1 < start2) {
                    LOG.info("{} in list 1 missing in list 2: {}-{}",featureKey, start1, end1);
                    i++;
                }
                // Fall: exon2 hat einen früheren Start, also fehlt dieses Exon vermutlich in list1
                else {
                    LOG.info("{} in list 2 missing in list 1: {}-{}", featureKey, start2, end2);
                    j++;
                }
            }

            while (i < sortedExons1.size()){
                GtfFeature exon1 = sortedExons1.get(i);
                LOG.info("Extra {} in list 1: {}-{}",featureKey, exon1.getBaseData().getStart(), exon1.getBaseData().getEnd());
                i++;
            }
            while (j < sortedExons2.size()){
                GtfFeature exon2 = sortedExons2.get(j);
                LOG.info("Extra {} in list 2: {}-{}",featureKey, exon2.getBaseData().getStart(), exon2.getBaseData().getEnd());
                j++;
            }

        }



    }
}
