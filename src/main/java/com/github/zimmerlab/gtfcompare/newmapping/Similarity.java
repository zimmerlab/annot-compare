package com.github.zimmerlab.gtfcompare.newmapping;

import com.github.kleinsamuel.gtfutils.GtfConfig;
import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.newmapping.model.CigarOp;
import com.github.zimmerlab.gtfcompare.newmapping.seqhomology.AlignmentUtil;
import com.github.zimmerlab.gtfcompare.newmapping.seqhomology.SeqHomologyUtil;
import com.github.zimmerlab.gtfcompare.utils.GenomeSequenceExtractor;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.IOException;
import java.util.*;

import static com.github.zimmerlab.gtfcompare.newmapping.MappingConstants.*;

public class Similarity {
    private final static Logger logger = LogManager.getLogger(Similarity.class);
    private static final SeqHomologyUtil seqHomologyUtil = new AlignmentUtil();

    public static List<CigarOp> structureCigar(TranscriptFeature feature, Set<String> allowedTypes, GenomeSequenceExtractor sequenceExtractor) {
        var result = new ArrayList<CigarOp>();
        var forwardStrand = feature.getBaseData().isForwardStrand();

        var comparator = Comparator.comparingInt((GtfFeature f) -> f.getBaseData().getStart());

        if (!forwardStrand) {
            comparator = comparator.reversed();
        }

        var features = feature.getFeatures().stream()
                .filter(f -> allowedTypes.contains(GtfConfig.getDefault(f.getBaseData().getType())))
                .sorted(comparator)
                .toList();
        
        for (var g : features) {
            var baseData = g.getBaseData();
            var type = GtfConfig.getDefault(baseData.getType());

            var start = baseData.getStart();
            var stop = baseData.getEnd();

            var len = stop - start + 1;
            String seq = null;

            if (features.size() == 1) {
                try {
                    seq = sequenceExtractor.getSequence(baseData.getContig(), start, stop);
                    if(!forwardStrand) seq = reverseComplement(seq);
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }

            result.add(new CigarOp(type, len, seq));
        }

        return result;
    }

    private static double cigarSimilarity(List<CigarOp> a, List<CigarOp> b) {
        if (a.isEmpty() && b.isEmpty()) return 0.0;

        if (a.isEmpty() || b.isEmpty()) return 0.0;

        var minSize = Math.min(a.size(), b.size());
        var maxSize = Math.max(a.size(), b.size());

        var similarityScore = 0.0;

        for (int i = 0; i < minSize; i++) {
            var opA = a.get(i);
            var opB = b.get(i);

            if (!opA.type().equals(opB.type())) continue;

            var lenA = opA.length();
            var lenB = opB.length();

            var maxLen = Math.max(lenA, lenB);

            if (maxLen == 0) {
                similarityScore += 1;
            } else {
                var relDiff = Math.abs(lenA - lenB) / (double) maxLen;
                var posSimilarity = 1.0 - relDiff;

                similarityScore += posSimilarity;
            }
        }

        return similarityScore / (double) maxSize;
    }

    private static double handleSingleFeature(List<CigarOp> a, List<CigarOp> b) {
        var opA = a.getFirst();
        var opB = b.getFirst();

        if (!opA.type().equals(opB.type())) return 0.0;

        var maxLength = Math.max(opA.length(), opB.length());

        if (maxLength == 0) return 0.0;

        var relDiff = Math.abs(opA.length() - opB.length()) / (double) maxLength;
        var baseSimilarity = 1.0 - relDiff;


        return baseSimilarity >= SIMILARITY_CUTOFF ? seqHomologyUtil.calculate(opA.seq(), opB.seq()) : 0;
    }

    public static boolean useSeqHomology = true;
    public static boolean isSimilar(List<CigarOp> a, List<CigarOp> b) {
        if (useSeqHomology && a.size() == 1 && b.size() == 1) {
            return handleSingleFeature(a, b) >= SEQ_HOMOLOGY_CUTOFF;
        }

        return cigarSimilarity(a, b) >= SIMILARITY_CUTOFF;
    }

    public static boolean isSimilarHomology(TranscriptFeature targetTranscript, TranscriptFeature queryTranscript, GenomeSequenceExtractor targetSequenceExtractor, GenomeSequenceExtractor querySequenceExtractor) {
        var targetSeq = new StringBuilder();
        var querySeq = new StringBuilder();


        var targetCds = targetTranscript.getFeatures().stream().filter(f -> GtfConfig.TYPE_CDS_SYNONYMS.contains(f.getBaseData().getType())).sorted(Comparator.comparingInt(f -> f.getBaseData().getStart())).toList();
        var queryCds = queryTranscript.getFeatures().stream().filter(f -> GtfConfig.TYPE_CDS_SYNONYMS.contains(f.getBaseData().getType())).sorted(Comparator.comparingInt(f -> f.getBaseData().getStart())).toList();

        boolean targetForward = targetCds.isEmpty() || targetCds.getFirst().getBaseData().isForwardStrand();
        boolean queryForward = queryCds.isEmpty() || queryCds.getFirst().getBaseData().isForwardStrand();

        for (var cds : targetCds) {
            var baseData = cds.getBaseData();

            var start = baseData.getStart();
            var stop = baseData.getEnd();
            var contig = baseData.getContig();

            String cdsSeq = null;
            try {
                cdsSeq = targetSequenceExtractor.getSequence(contig, start, stop);

            } catch (IOException e) {
                logger.warn(e.getMessage());
                continue;
            }

            if (cdsSeq == null) continue;
            targetSeq.append(cdsSeq);
        }

        for (var cds : queryCds) {
            var baseData = cds.getBaseData();

            var start = baseData.getStart();
            var stop = baseData.getEnd();
            var contig = baseData.getContig();

            String cdsSeq = null;
            try {
                cdsSeq = querySequenceExtractor.getSequence(contig, start, stop);
            } catch (IOException e) {
                logger.warn(e.getMessage());
                continue;
            }

            if (cdsSeq == null) continue;
            querySeq.append(cdsSeq);
        }

        var targetSeqString = targetForward ? targetSeq.toString() : reverseComplement(targetSeq.toString());
        var querySeqString = queryForward ? querySeq.toString() : reverseComplement(querySeq.toString());
        return targetSeqString.equals(querySeqString);
    }

    private static String reverseComplement(String seq) {
        StringBuilder revComp = new StringBuilder(seq.length());

        for (int i = seq.length() - 1; i >= 0; i--) {
            char base = seq.charAt(i);
            switch (base) {
                case 'A':
                    revComp.append('T');
                    break;
                case 'T':
                    revComp.append('A');
                    break;
                case 'C':
                    revComp.append('G');
                    break;
                case 'G':
                    revComp.append('C');
                    break;
                case 'a':
                    revComp.append('t');
                    break;
                case 't':
                    revComp.append('a');
                    break;
                case 'c':
                    revComp.append('g');
                    break;
                case 'g':
                    revComp.append('c');
                    break;
                default:
                    revComp.append('N');
            }
        }

        return revComp.toString();
    }
}
