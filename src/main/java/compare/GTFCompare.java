package compare;

import com.github.kleinsamuel.gtfutils.feature.GeneFeature;
import com.github.kleinsamuel.gtfutils.feature.GtfBaseData;
import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.GenomeSequenceExtractor;
import com.github.zimmerlab.gtfcompare.model.TranscriptPair;
import com.github.zimmerlab.gtfcompare.model.comparison.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutionException;

// TODO analyse und auslagern
public class GTFCompare {
    public static final Logger LOG = LoggerFactory.getLogger(GTFCompare.class);

    private static GenomeSequenceExtractor sequenceExtractor1;
    private static GenomeSequenceExtractor sequenceExtractor2;

    public static ComparisonResult compare(String geneId1, String geneId2,
                                           GeneFeature gene1, GeneFeature gene2,
                                           GenomeSequenceExtractor genomeSequenceExtractor1,
                                           GenomeSequenceExtractor genomeSequenceExtractor2) {

        sequenceExtractor1 = genomeSequenceExtractor1;
        sequenceExtractor2 = genomeSequenceExtractor2;

        ComparisonResult comparisonResult = new ComparisonResult();
        comparisonResult.setGeneId1(geneId1);
        comparisonResult.setGeneId2(geneId2);

        compareGenes(gene1, gene2, comparisonResult);

        if(!comparisonResult.areSameGene())
            return comparisonResult;

        compareTranscripts(gene1, gene2, comparisonResult, comparisonResult.getGeneComparison().getSequenceComparison().isSameSequence());

        return comparisonResult;
    }

    private static void compareGenes(GeneFeature gene1, GeneFeature gene2, ComparisonResult comparisonResult){
        var baseData1 = gene1.getBaseData();
        var baseData2 = gene2.getBaseData();
        var geneComparisonResult = comparisonResult.getGeneComparison();


        var start1 = baseData1.getStart();
        var stop1 = baseData1.getEnd();

        var start2 = baseData2.getStart();
        var stop2 = baseData2.getEnd();

        if((stop1 - start1) != (stop2 - start2)){
            geneComparisonResult.setDifferentLength(true);
            comparisonResult.setAreSameGene(false);
            return;
        }

        if (baseData1.isForwardStrand() != baseData2.isForwardStrand()) {
            geneComparisonResult.setStrandDifferent(true);
            comparisonResult.setAreSameGene(false);
            geneComparisonResult.addMessage("Strand changed from " + baseData1.isForwardStrand() + " to " + baseData2.isForwardStrand());
            return;
        }

        if (start1 != start2) {
            geneComparisonResult.setStartDifferent(true);
            geneComparisonResult.addMessage("Start position not the same");
        }

        if (stop1 != stop2) {
            geneComparisonResult.setStopDifferent(true);
            geneComparisonResult.addMessage("Stop position not the same");
        }

        try {
            String seq1 = sequenceExtractor1.getSequence(baseData1.getContig(), baseData1.getStart(), baseData1.getEnd());
            String seq2 = sequenceExtractor2.getSequence(baseData2.getContig(), baseData2.getStart(), baseData2.getEnd());


            if (!seq1.equals(seq2)) {
                comparisonResult.setAreSameGene(false);
                var sequenceComparison = geneComparisonResult.getSequenceComparison();
                sequenceComparison.setSameSequence(false);
                sequenceComparison.setSeq1(seq1);
                sequenceComparison.setSeq2(seq2);
                geneComparisonResult.addMessage("Gene sequence not the same");
                if (seq1.length() != seq2.length()) {
                    geneComparisonResult.addMessage("Different sequence lengths: seq1 = " + seq1.length() + ", seq2 = " + seq2.length());
                }
                return;
                // TODO alignment?
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private static List<TranscriptPair> getTranscriptPairs(GeneFeature g1, GeneFeature g2, ComparisonResult comparisonResult) {
        var transcriptsMap1 = new HashMap<String, TranscriptFeature>();
        var transcriptsMap2 = new HashMap<String, TranscriptFeature>();

        for (TranscriptFeature transcript : g1.getTranscripts()) {
            transcriptsMap1.put(transcript.getTranscriptId(), transcript);
        }

        for (TranscriptFeature transcript : g2.getTranscripts()) {
            transcriptsMap2.put(transcript.getTranscriptId(), transcript);
        }

        var allTranscriptIds = new HashSet<String>();
        allTranscriptIds.addAll(transcriptsMap1.keySet());
        allTranscriptIds.addAll(transcriptsMap2.keySet());

        var transcriptPairs = new ArrayList<TranscriptPair>();
        for (var transcriptId : allTranscriptIds) {
            var transcriptComparisonResult = new TranscriptComparisonResult();
            TranscriptFeature t1 = transcriptsMap1.get(transcriptId);
            TranscriptFeature t2 = transcriptsMap2.get(transcriptId);

            if (t1 == null) {
                comparisonResult.setAreSameGene(false);
                transcriptComparisonResult.setTranscriptMissingInGene1(true);
                transcriptComparisonResult.setTranscriptId(transcriptId);
                comparisonResult.addTranscriptComparison(transcriptComparisonResult);
                return null;
            } else if (t2 == null) {
                comparisonResult.setAreSameGene(false);
                transcriptComparisonResult.setTranscriptMissingInGene2(true);
                transcriptComparisonResult.setTranscriptId(transcriptId);
                comparisonResult.addTranscriptComparison(transcriptComparisonResult);
                return null;
            }

            transcriptPairs.add(new TranscriptPair(t1, t2, transcriptComparisonResult));
        }

        return transcriptPairs;
    }

    private static boolean compareTranscriptDetails(TranscriptFeature t1,
                                                    TranscriptFeature t2,
                                                    boolean isSameGeneSequence,
                                                    TranscriptComparisonResult transcriptComparisonResult, ComparisonResult comparisonResult) {
        GtfBaseData baseData1 = t1.getBaseData();
        GtfBaseData baseData2 = t2.getBaseData();

        var start1 = baseData1.getStart();
        var stop1 = baseData1.getEnd();

        var start2 = baseData2.getStart();
        var stop2 = baseData2.getEnd();

        boolean isStartDifferent = start1 != start2;
        boolean isStopDifferent = stop1 != stop2;

        if (isStartDifferent) {
            transcriptComparisonResult.setStartDifferent(true);
        }
        if (isStopDifferent) {
            transcriptComparisonResult.setStopDifferent(true);
        }

        if((stop1 - start1) != (stop2 - start2)){
            transcriptComparisonResult.setLengthDifferent(true);
            return false;
        }

        boolean isSameTranscriptSequence = true;

        if (!isSameGeneSequence || isStartDifferent || isStopDifferent) {
            var sequenceComparisonResult = isSameSequence(
                    baseData1.getContig(),
                    baseData1.getStart(),
                    baseData1.getEnd(),
                    baseData2.getStart(),
                    baseData2.getEnd());
            if (sequenceComparisonResult.getSeq1() != null) {
                LOG.info("Transcript sequence not the same");
                isSameTranscriptSequence = false;
                comparisonResult.setAreSameGene(false);
                return false;
            }
        }

        return isSameTranscriptSequence;
    }

    private static void compareTranscripts(GeneFeature g1, GeneFeature g2, ComparisonResult comparisonResult, boolean isSameGeneSequence) {
        var transcriptPairs = getTranscriptPairs(g1, g2, comparisonResult);
        if(!comparisonResult.areSameGene()){
            return;
        }
        for (var transcripts : transcriptPairs) {
            TranscriptFeature t1 = transcripts.getTranscript1();
            TranscriptFeature t2 = transcripts.getTranscript2();
            var transcriptComparisonResult = transcripts.getTranscriptComparisonResult();

            try {
                transcriptComparisonResult.setTranscriptId(t1.getTranscriptId());
            } catch (Exception e) {
                try{
                    transcriptComparisonResult.setTranscriptId(t2.getTranscriptId());
                } catch (Exception ignored){}
            }

            comparisonResult.addTranscriptComparison(transcriptComparisonResult);
            boolean isSameTranscriptSequence = false;
            if (t1 != null && t2 != null) {
                isSameTranscriptSequence = compareTranscriptDetails(t1, t2, isSameGeneSequence, transcriptComparisonResult, comparisonResult);
                if(!comparisonResult.areSameGene())
                    return;
                compareFeatures(t1, t2, transcriptComparisonResult, isSameTranscriptSequence, comparisonResult);
                if(!comparisonResult.areSameGene())
                    return;
            }
        }
    }

    private static void compareFeatures(TranscriptFeature t1, TranscriptFeature t2, TranscriptComparisonResult transcriptComparisonResult, boolean isSameTranscriptSequence, ComparisonResult comparisonResult) {
        var featureMap1 = new HashMap<String, List<GtfFeature>>();
        var featureMap2 = new HashMap<String, List<GtfFeature>>();

        for (var feature : t1.getFeatures()) {
            featureMap1.computeIfAbsent(feature.getBaseData().getType(), key -> new ArrayList<>()).add(feature);
        }
        for (var feature : t2.getFeatures()) {
            featureMap2.computeIfAbsent(feature.getBaseData().getType(), key -> new ArrayList<>()).add(feature);
        }
        compareFeaturePositions(t1, t2, featureMap1, featureMap2, transcriptComparisonResult, isSameTranscriptSequence, comparisonResult);
    }

    private static void compareFeaturePositions(TranscriptFeature t1, TranscriptFeature t2,
                                                   Map<String, List<GtfFeature>> featureMap1,
                                                   Map<String, List<GtfFeature>> featureMap2,
                                                   TranscriptComparisonResult transcriptComparisonResult,
                                                   boolean isSameTranscriptSequence, ComparisonResult comparisonResult) {
        var featureComparisons = transcriptComparisonResult.getFeatureComparisons();
        // TODO length -> seq -> position
        for (var featureKey : featureMap2.keySet()) {
            var featureComparisonResult = new FeatureComparisonResult();
            featureComparisonResult.setFeatureType(featureKey);

            if (!featureMap1.containsKey(featureKey)) {
                featureComparisonResult.setMissingInTranscript1(true);
                String msg = featureKey + " not in " + t1.getTranscriptId();
                featureComparisons.add(featureComparisonResult);
                LOG.info(msg);
                comparisonResult.setAreSameGene(false);
                return;
            }
        }

        for (var featureKey : featureMap1.keySet()) {
            var featureComparisonResult = new FeatureComparisonResult();
            featureComparisonResult.setFeatureType(featureKey);

            // If the feature type is missing from the second transcript, record it.
            if (!featureMap2.containsKey(featureKey)) {
                featureComparisonResult.setMissingInTranscript2(true);
                String msg = featureKey + " not in " + t2.getTranscriptId();
                //transcriptComparisonResult.addMessage(msg);
                //featureComparisonResult.addMessage(msg);
                featureComparisons.add(featureComparisonResult);
                LOG.info(msg);
                comparisonResult.setAreSameGene(false);
                return;
            }

            List<GtfFeature> sortedFeatures1 = new ArrayList<>(featureMap1.get(featureKey));
            List<GtfFeature> sortedFeatures2 = new ArrayList<>(featureMap2.get(featureKey));

            sortedFeatures1.sort(
                    Comparator.comparingInt((GtfFeature f) -> f.getBaseData().getStart())
                            .thenComparingInt(f -> f.getBaseData().getEnd())
            );
            sortedFeatures2.sort(
                    Comparator.comparingInt((GtfFeature f) -> f.getBaseData().getStart())
                            .thenComparingInt(f -> f.getBaseData().getEnd())
            );

            int i = 0, j = 0;
            while (i < sortedFeatures1.size() && j < sortedFeatures2.size()) {
                GtfFeature feature1 = sortedFeatures1.get(i);
                GtfFeature feature2 = sortedFeatures2.get(j);

                var baseData1 = feature1.getBaseData();
                var baseData2 = feature2.getBaseData();

                var currentContig = baseData1.getContig();

                int start1 = baseData1.getStart();
                int start2 = baseData2.getStart();
                int end1 = baseData1.getEnd();
                int end2 = baseData2.getEnd();

                var length1 = end1 - start1;
                var length2 = end2 - start2;
                var regionComparison = new RegionComparison(start1, end1, start2, end2);
                featureComparisonResult.addRegionComparison(regionComparison);
                featureComparisons.add(featureComparisonResult);

                i++;
                j++;
                if(length1 != length2){
                    regionComparison.setLengthDifferenceFound(true);
                    comparisonResult.setAreSameGene(false);
                    return;
                }

                if(start1 != start2 || end1 != end2){
                    regionComparison.setPositionDifferenceFound(true);
                }

                if(!isSameTranscriptSequence || (regionComparison.isPositionDifferenceFound() || regionComparison.isLengthDifferenceFound())){
                    try{
                        var seq1 = sequenceExtractor1.getSequence(currentContig, start1, end1);
                        var seq2 = sequenceExtractor2.getSequence(currentContig, start2, end2);

                        if(!seq1.equals(seq2)){
                            regionComparison.setSequenceDifferenceFound(true);
                            comparisonResult.setAreSameGene(false);
                            return;
                        }

                    } catch (Exception ignored){}
                }
            }

            /*// Process any leftover features in sortedFeatures1.
            while (i < sortedFeatures1.size()) {
                GtfFeature feature1 = sortedFeatures1.get(i);
                var start = feature1.getBaseData().getStart();
                var stop = feature1.getBaseData().getEnd();
                String msg = featureKey + " in " + t2.getTranscriptId() +
                        " missing in list 1: " + start + "-" + stop;
                LOG.info(msg);
                featureComparisonResult.addRegionComparison(new RegionComparison(start, stop, -1, -1, true));
                //featureComparisonResult.addMessage(msg);
                i++;
            }

            // Process any leftover features in sortedFeatures2.
            while (j < sortedFeatures2.size()) {
                GtfFeature feature2 = sortedFeatures2.get(j);
                var start = feature2.getBaseData().getStart();
                var stop = feature2.getBaseData().getEnd();
                String msg = featureKey + " in " + t2.getTranscriptId() +
                        " missing in list 2: " + start + "-" + stop;
                LOG.info(msg);
                //featureComparisonResult.setMissingInTranscript1(true);
                featureComparisonResult.addRegionComparison(new RegionComparison(-1, -1, start, stop, true));
                j++;
            }*/

        }
    }

    private static void compareFeaturePositionsOld(TranscriptFeature t1, TranscriptFeature t2,
                                                Map<String, List<GtfFeature>> featureMap1,
                                                Map<String, List<GtfFeature>> featureMap2,
                                                TranscriptComparisonResult transcriptComparisonResult,
                                                boolean isSameTranscriptSequence) {
        var featureComparisons = transcriptComparisonResult.getFeatureComparisons();
        // TODO length -> seq -> position
        for (var featureKey : featureMap2.keySet()) {
            var featureComparisonResult = new FeatureComparisonResult();
            featureComparisonResult.setFeatureType(featureKey);

            if (!featureMap1.containsKey(featureKey)) {
                featureComparisonResult.setMissingInTranscript1(true);
                String msg = featureKey + " not in " + t1.getTranscriptId();
                featureComparisons.add(featureComparisonResult);
                LOG.info(msg);
            }
        }

        for (var featureKey : featureMap1.keySet()) {
            var featureComparisonResult = new FeatureComparisonResult();
            featureComparisonResult.setFeatureType(featureKey);

            // If the feature type is missing from the second transcript, record it.
            if (!featureMap2.containsKey(featureKey)) {
                featureComparisonResult.setMissingInTranscript2(true);
                String msg = featureKey + " not in " + t2.getTranscriptId();
                //transcriptComparisonResult.addMessage(msg);
                //featureComparisonResult.addMessage(msg);
                featureComparisons.add(featureComparisonResult);
                LOG.info(msg);
                continue;
            }

            List<GtfFeature> sortedFeatures1 = new ArrayList<>(featureMap1.get(featureKey));
            List<GtfFeature> sortedFeatures2 = new ArrayList<>(featureMap2.get(featureKey));

            sortedFeatures1.sort(
                    Comparator.comparingInt((GtfFeature f) -> f.getBaseData().getStart())
                            .thenComparingInt(f -> f.getBaseData().getEnd())
            );
            sortedFeatures2.sort(
                    Comparator.comparingInt((GtfFeature f) -> f.getBaseData().getStart())
                            .thenComparingInt(f -> f.getBaseData().getEnd())
            );

            int i = 0, j = 0;
            while (i < sortedFeatures1.size() && j < sortedFeatures2.size()) {
                GtfFeature feature1 = sortedFeatures1.get(i);
                GtfFeature feature2 = sortedFeatures2.get(j);

                var baseData1 = feature1.getBaseData();
                var baseData2 = feature2.getBaseData();

                int start1 = baseData1.getStart();
                int start2 = baseData2.getStart();
                int end1 = baseData1.getEnd();
                int end2 = baseData2.getEnd();
                var regionComparison = new RegionComparison(start1, end1, start2, end2);
                // Case 1: Both start and end are identical.
                if (start1 == start2 && end1 == end2) {
                    i++;
                    j++;

                    if (!isSameTranscriptSequence) {
                        var sequenceComparisonResult = isSameSequence(baseData1.getContig(), start1, end1, start2, end2);
                        if (sequenceComparisonResult.getSeq1() != null) {
                            LOG.info("{} has not the same sequence start: {}, stop: {}", featureKey, start1, end1);
                            regionComparison.setSequenceDifferenceFound(true);
                        }
                    }

                }
                // Case 2: Same start but different end positions.
                else if (start1 == start2) {
                    // TODO evtl. einfach als missing markieren?
                    String msg = featureKey + " end changed: list 1: " + start1 + "-" + end1 +
                            ", list 2: " + start2 + "-" + end2;
                    LOG.info(msg);
                    regionComparison.setPositionDifferenceFound(true);
                    //featureComparisonResult.addMessage(msg);
                    i++;
                    j++;

                    var sequenceComparisonResult = isSameSequence(baseData1.getContig(), start1, end1, start2, end2);
                    if (sequenceComparisonResult.getSeq1() != null) {
                        LOG.info("{} has not the same sequence start: {}, stop: {}", featureKey, start1, end1);
                        regionComparison.setSequenceDifferenceFound(true);
                    }
                }
                // Case 3: Same end but different start positions.
                else if (end1 == end2) {
                    // TODO evtl. einfach als missing markieren?
                    String msg = featureKey + " start changed: list 1: " + start1 + "-" + end1 +
                            ", list 2: " + start2 + "-" + end2;
                    LOG.info(msg);
                    regionComparison.setPositionDifferenceFound(true);
                    //featureComparisonResult.addRegionComparison(regionComparison);
                    //featureComparisonResult.addMessage(msg);
                    i++;
                    j++;

                    var sequenceComparisonResult = isSameSequence(baseData1.getContig(), start1, end1, start2, end2);
                    if (sequenceComparisonResult.getSeq1() != null) {
                        LOG.info("{} has not the same sequence start: {}, stop: {}", featureKey, start1, end1);
                        regionComparison.setSequenceDifferenceFound(true);
                    }
                }
                // Case 4: Feature from list 1 starts earlier, indicating it is missing in t2.
                else if (start1 < start2) {
                    String msg = featureKey + " in " + t1.getTranscriptId() +
                            " missing in list 2: " + start1 + "-" + end1;
                    LOG.info(msg);
                    regionComparison.setStart2(-1);
                    regionComparison.setEnd2(-1);
                    regionComparison.setPositionDifferenceFound(true);
                    //featureComparisonResult.addMessage(msg);
                    i++;
                }
                // Case 5: Feature from list 2 starts earlier, indicating it is missing in t1.
                else {
                    String msg = featureKey + " in " + t2.getTranscriptId() +
                            " missing in list 1: " + start2 + "-" + end2;
                    LOG.info(msg);
                    regionComparison.setStart1(-1);
                    regionComparison.setEnd1(-1);
                    regionComparison.setPositionDifferenceFound(true);
                    //featureComparisonResult.addMessage(msg);
                    j++;
                }

                featureComparisonResult.addRegionComparison(regionComparison);
            }

            // Process any leftover features in sortedFeatures1.
            while (i < sortedFeatures1.size()) {
                GtfFeature feature1 = sortedFeatures1.get(i);
                var start = feature1.getBaseData().getStart();
                var stop = feature1.getBaseData().getEnd();
                String msg = featureKey + " in " + t2.getTranscriptId() +
                        " missing in list 1: " + start + "-" + stop;
                LOG.info(msg);
                featureComparisonResult.addRegionComparison(new RegionComparison(start, stop, -1, -1, true));
                //featureComparisonResult.addMessage(msg);
                i++;
            }

            // Process any leftover features in sortedFeatures2.
            while (j < sortedFeatures2.size()) {
                GtfFeature feature2 = sortedFeatures2.get(j);
                var start = feature2.getBaseData().getStart();
                var stop = feature2.getBaseData().getEnd();
                String msg = featureKey + " in " + t2.getTranscriptId() +
                        " missing in list 2: " + start + "-" + stop;
                LOG.info(msg);
                //featureComparisonResult.setMissingInTranscript1(true);
                featureComparisonResult.addRegionComparison(new RegionComparison(-1, -1, start, stop, true));
                j++;
            }

            featureComparisons.add(featureComparisonResult);
        }
    }

    private static SequenceComparisonResult isSameSequence(String chr, int start1, int stop1, int start2, int stop2) {
        try {
            String seq1 = sequenceExtractor1.getSequence(chr, start1, stop1);
            String seq2 = sequenceExtractor2.getSequence(chr, start2, stop2);
            var sequenceComparisonResult = new SequenceComparisonResult();

            if (!seq1.equals(seq2)) {
                sequenceComparisonResult.setSeq1(seq1);
                sequenceComparisonResult.setSeq2(seq2);
                // TODO add differences
                // TODO add type of difference(s)
                // TODO alignment?
            }

            return sequenceComparisonResult;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

}