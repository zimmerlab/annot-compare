package compare;

import com.github.kleinsamuel.gtfutils.feature.GeneFeature;
import com.github.kleinsamuel.gtfutils.feature.GtfBaseData;
import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.utils.Constants;
import com.github.zimmerlab.gtfcompare.utils.GenomeSequenceExtractor;
import com.github.zimmerlab.gtfcompare.model.TranscriptPair;
import com.github.zimmerlab.gtfcompare.model.comparison.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.util.StopWatch;

import java.io.IOException;
import java.util.*;

public class GTFCompare {
    public static final Logger LOG = LoggerFactory.getLogger(GTFCompare.class);

    private static GenomeSequenceExtractor sequenceExtractor1;
    private static GenomeSequenceExtractor sequenceExtractor2;
    private static final List<StopWatch> stopWatches = Constants.STOP_WATCHES;

    // TODO do I need to check the sequence if the position changed but length and sequence on the lvl above stayed the same?
    public static ComparisonResult compare(String targetGeneId, String queryGeneId,
                                           GeneFeature targetGene, GeneFeature queryGene,
                                           GenomeSequenceExtractor targetGenomeSequenceExtractor,
                                           GenomeSequenceExtractor queryGenomeSequenceExtractor) {

        sequenceExtractor1 = targetGenomeSequenceExtractor;
        sequenceExtractor2 = queryGenomeSequenceExtractor;

        ComparisonResult comparisonResult = new ComparisonResult();
        comparisonResult.setTargetGeneId(targetGeneId);
        comparisonResult.setQueryGeneId(queryGeneId);

        compareGenes(targetGene, queryGene, comparisonResult);

        compareTranscripts(targetGene, queryGene, comparisonResult, comparisonResult.getGeneComparison().getSequenceComparison().isSameSequence());

        return comparisonResult;
    }

    private static void compareGenes(GeneFeature targetGene, GeneFeature queryGene, ComparisonResult comparisonResult) {
        var stopWatch = new StopWatch();
        stopWatches.add(stopWatch);
        stopWatch.start("compareGene");

        var targetBaseData = targetGene.getBaseData();
        var queryBaseData = queryGene.getBaseData();
        var geneComparisonResult = comparisonResult.getGeneComparison();


        var start1 = targetBaseData.getStart();
        var stop1 = targetBaseData.getEnd();

        var start2 = queryBaseData.getStart();
        var stop2 = queryBaseData.getEnd();

        if ((stop1 - start1) != (stop2 - start2)) {
            geneComparisonResult.setLengthDifferent(true);
            comparisonResult.setAreSameGene(false);
        }

        if (targetBaseData.isForwardStrand() != queryBaseData.isForwardStrand()) {
            geneComparisonResult.setStrandDifferent(true);
            comparisonResult.setAreSameGene(false);
            geneComparisonResult.addMessage("Strand changed from " + targetBaseData.isForwardStrand() + " to " + queryBaseData.isForwardStrand());
        }

        if (start1 != start2) {
            geneComparisonResult.setStartDifferent(true);
            geneComparisonResult.addMessage("Start position not the same");
        }

        if (stop1 != stop2) {
            geneComparisonResult.setStopDifferent(true);
            geneComparisonResult.addMessage("Stop position not the same");
        }

        if (geneComparisonResult.isDifferentLength()) {
            comparisonResult.setAreSameGene(false);
            var sequenceComparison = geneComparisonResult.getSequenceComparison();
            sequenceComparison.setSameSequence(false);
            geneComparisonResult.addMessage("Gene sequence not the same");
            return;
        }

        stopWatch.stop();
        stopWatch.start("compareGene_SequenceExtraction");
        try {
            String seq1 = sequenceExtractor1.getSequence(targetBaseData.getContig(), targetBaseData.getStart(), targetBaseData.getEnd());
            String seq2 = sequenceExtractor2.getSequence(queryBaseData.getContig(), queryBaseData.getStart(), queryBaseData.getEnd());

            if (!seq1.equals(seq2)) {
                comparisonResult.setAreSameGene(false);
                var sequenceComparison = geneComparisonResult.getSequenceComparison();
                sequenceComparison.setSameSequence(false);
                geneComparisonResult.addMessage("Gene sequence not the same");
                // TODO alignment?
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        stopWatch.stop();
    }

    private static List<TranscriptPair> getTranscriptPairs(GeneFeature targetGene, GeneFeature queryGene, ComparisonResult comparisonResult) {
        var stopWatch = new StopWatch();
        stopWatches.add(stopWatch);
        stopWatch.start("getTranscriptPairs");
        var transcriptsMap1 = new HashMap<String, TranscriptFeature>();
        var transcriptsMap2 = new HashMap<String, TranscriptFeature>();

        for (TranscriptFeature transcript : targetGene.getTranscripts()) {
            transcriptsMap1.put(transcript.getTranscriptId(), transcript);
        }

        for (TranscriptFeature transcript : queryGene.getTranscripts()) {
            transcriptsMap2.put(transcript.getTranscriptId(), transcript);
        }

        var allTranscriptIds = new HashSet<String>();
        allTranscriptIds.addAll(transcriptsMap1.keySet());
        allTranscriptIds.addAll(transcriptsMap2.keySet());

        var transcriptPairs = new ArrayList<TranscriptPair>();
        for (var transcriptId : allTranscriptIds) {
            var transcriptComparisonResult = new TranscriptComparisonResult();
            comparisonResult.addTranscriptComparison(transcriptComparisonResult);
            TranscriptFeature t1 = transcriptsMap1.get(transcriptId);
            TranscriptFeature t2 = transcriptsMap2.get(transcriptId);

            if (t1 == null) {
                comparisonResult.setAreSameGene(false);
                transcriptComparisonResult.setTranscriptMissingInTargetGene(true);
                transcriptComparisonResult.setTargetTranscriptId(transcriptId);
            } else if (t2 == null) {
                comparisonResult.setAreSameGene(false);
                transcriptComparisonResult.setTranscriptMissingInQueryGene(true);
                transcriptComparisonResult.setQueryTranscriptId(transcriptId);
            }

            transcriptPairs.add(new TranscriptPair(t1, t2, transcriptComparisonResult));
        }

        stopWatch.stop();
        return transcriptPairs;
    }

    private static boolean compareTranscriptDetails(TranscriptFeature t1,
                                                    TranscriptFeature t2,
                                                    boolean isSameGeneSequence,
                                                    TranscriptComparisonResult transcriptComparisonResult, ComparisonResult comparisonResult) {
        var stopWatch = new StopWatch();
        stopWatches.add(stopWatch);
        stopWatch.start("compareTranscriptDetails");
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

        if ((stop1 - start1) != (stop2 - start2)) {
            transcriptComparisonResult.setLengthDifferent(true);
            transcriptComparisonResult.setSequenceDifferent(true);
            comparisonResult.setAreSameGene(false);
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
            if (sequenceComparisonResult.getTargetSeq() != null) {
                LOG.info("Transcript sequence not the same");
                isSameTranscriptSequence = false;
                comparisonResult.setAreSameGene(false);
                transcriptComparisonResult.setSequenceDifferent(true);
            }
        }
        stopWatch.stop();
        return isSameTranscriptSequence;
    }

    private static void compareTranscripts(GeneFeature g1, GeneFeature g2, ComparisonResult comparisonResult, boolean isSameGeneSequence) {
        var transcriptPairs = getTranscriptPairs(g1, g2, comparisonResult);

        for (var transcripts : transcriptPairs) {
            TranscriptFeature t1 = transcripts.getTarget();
            TranscriptFeature t2 = transcripts.getQuery();
            var transcriptComparisonResult = transcripts.getTranscriptComparisonResult();

            try {
                transcriptComparisonResult.setQueryTranscriptId(t1.getTranscriptId());
            } catch (Exception e) {
                try {
                    transcriptComparisonResult.setTargetTranscriptId(t2.getTranscriptId());
                } catch (Exception ignored) {
                }
            }

            boolean isSameTranscriptSequence = false;
            if (t1 != null && t2 != null) {
                isSameTranscriptSequence = compareTranscriptDetails(t1, t2, isSameGeneSequence, transcriptComparisonResult, comparisonResult);
                compareFeatures(t1, t2, transcriptComparisonResult, isSameTranscriptSequence, comparisonResult);
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

        var stopWatch = new StopWatch();
        stopWatches.add(stopWatch);
        stopWatch.start("compareFeatures");
        var featureComparisons = transcriptComparisonResult.getFeatureComparisons();
        for (var featureKey : featureMap2.keySet()) {
            var featureComparisonResult = new FeatureComparisonResult();
            featureComparisonResult.setFeatureType(featureKey);

            if (!featureMap1.containsKey(featureKey)) {
                featureComparisonResult.setMissingInTargetTranscript(true);
                String msg = featureKey + " not in " + t1.getTranscriptId();
                featureComparisons.add(featureComparisonResult);
                LOG.info(msg);
                comparisonResult.setAreSameGene(false);
            }
        }

        for (var featureKey : featureMap1.keySet()) {
            var featureComparisonResult = new FeatureComparisonResult();
            featureComparisonResult.setFeatureType(featureKey);

            // If the feature type is missing from the second transcript, record it.
            if (!featureMap2.containsKey(featureKey)) {
                featureComparisonResult.setMissingInQueryTranscript(true);
                String msg = featureKey + " not in " + t2.getTranscriptId();
                //transcriptComparisonResult.addMessage(msg);
                //featureComparisonResult.addMessage(msg);
                featureComparisons.add(featureComparisonResult);
                LOG.info(msg);
                comparisonResult.setAreSameGene(false);
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

                i++;
                j++;

                var baseData1 = feature1.getBaseData();
                var baseData2 = feature2.getBaseData();

                var currentContig = baseData1.getContig();

                int start1 = baseData1.getStart();
                int start2 = baseData2.getStart();
                int end1 = baseData1.getEnd();
                int end2 = baseData2.getEnd();

                var length1 = end1 - start1;
                var length2 = end2 - start2;
                var regionComparison = new RegionComparisonResult(start1, end1, start2, end2);
                featureComparisonResult.addRegionComparison(regionComparison);
                featureComparisons.add(featureComparisonResult);

                if (length1 != length2) {
                    regionComparison.setLengthDifferent(true);
                    regionComparison.setSequenceDifferenceFound(true);
                    comparisonResult.setAreSameGene(false);
                }

                if (start1 != start2 || end1 != end2) {
                    regionComparison.setPositionDifferenceFound(true);
                }

                if (!regionComparison.isSequenceDifferenceFound() && (!isSameTranscriptSequence || regionComparison.isLengthDifferenceFound())) {
                    try {
                        stopWatch.stop();
                        stopWatch.start("compareFeatures_Sequence");
                        var seq1 = sequenceExtractor1.getSequence(currentContig, start1, end1);
                        var seq2 = sequenceExtractor2.getSequence(currentContig, start2, end2);

                        if (!seq1.equals(seq2)) {
                            regionComparison.setSequenceDifferenceFound(true);
                            comparisonResult.setAreSameGene(false);
                        }

                        stopWatch.stop();
                        stopWatch.start("compareFeatures");

                    } catch (Exception ignored) {
                        if (stopWatch.isRunning()){
                            stopWatch.stop();
                            stopWatch.start("compareFeatures");
                        }
                    }
                }
            }

            // Process any leftover features in sortedFeatures1.
            while (i < sortedFeatures1.size()) {
                GtfFeature feature1 = sortedFeatures1.get(i);
                var start = feature1.getBaseData().getStart();
                var stop = feature1.getBaseData().getEnd();
                String msg = featureKey + " in " + t2.getTranscriptId() +
                        " missing in list 1: " + start + "-" + stop;
                LOG.info(msg);
                var regionComparison = new RegionComparisonResult(start, stop, -1, -1);
                regionComparison.setMissingInQueryFile(true);
                featureComparisonResult.addRegionComparison(regionComparison);
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
                var regionComparison = new RegionComparisonResult(-1, -1, start, stop);
                regionComparison.setMissingInTargetFile(true);
                featureComparisonResult.addRegionComparison(regionComparison);
                j++;
            }
            stopWatch.stop();
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
                featureComparisonResult.setMissingInTargetTranscript(true);
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
                featureComparisonResult.setMissingInQueryTranscript(true);
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
                var regionComparison = new RegionComparisonResult(start1, end1, start2, end2);
                // Case 1: Both start and end are identical.
                if (start1 == start2 && end1 == end2) {
                    i++;
                    j++;

                    if (!isSameTranscriptSequence) {
                        var sequenceComparisonResult = isSameSequence(baseData1.getContig(), start1, end1, start2, end2);
                        if (sequenceComparisonResult.getTargetSeq() != null) {
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
                    if (sequenceComparisonResult.getTargetSeq() != null) {
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
                    if (sequenceComparisonResult.getTargetSeq() != null) {
                        LOG.info("{} has not the same sequence start: {}, stop: {}", featureKey, start1, end1);
                        regionComparison.setSequenceDifferenceFound(true);
                    }
                }
                // Case 4: Feature from list 1 starts earlier, indicating it is missing in t2.
                else if (start1 < start2) {
                    String msg = featureKey + " in " + t1.getTranscriptId() +
                            " missing in list 2: " + start1 + "-" + end1;
                    LOG.info(msg);
                    regionComparison.setQueryStart(-1);
                    regionComparison.setQueryEnd(-1);
                    regionComparison.setPositionDifferenceFound(true);
                    //featureComparisonResult.addMessage(msg);
                    i++;
                }
                // Case 5: Feature from list 2 starts earlier, indicating it is missing in t1.
                else {
                    String msg = featureKey + " in " + t2.getTranscriptId() +
                            " missing in list 1: " + start2 + "-" + end2;
                    LOG.info(msg);
                    regionComparison.setTargetStart(-1);
                    regionComparison.setTargetEnd(-1);
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
                featureComparisonResult.addRegionComparison(new RegionComparisonResult(start, stop, -1, -1, true));
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
                featureComparisonResult.addRegionComparison(new RegionComparisonResult(-1, -1, start, stop, true));
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
                sequenceComparisonResult.setTargetSeq(seq1);
                sequenceComparisonResult.setQuerySeq(seq2);
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