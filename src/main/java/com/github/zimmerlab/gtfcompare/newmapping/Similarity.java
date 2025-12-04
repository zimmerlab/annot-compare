package com.github.zimmerlab.gtfcompare.newmapping;

import com.github.kleinsamuel.gtfutils.GtfConfig;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.newmapping.model.CigarOp;
import com.github.zimmerlab.gtfcompare.newmapping.seqhomology.AlignmentUtil;
import com.github.zimmerlab.gtfcompare.newmapping.seqhomology.SeqHomologyUtil;
import com.github.zimmerlab.gtfcompare.utils.GenomeSequenceExtractor;

import java.io.IOException;
import java.util.*;

import static com.github.zimmerlab.gtfcompare.newmapping.MappingConstants.*;

public class Similarity {
    private static final SeqHomologyUtil seqHomologyUtil = new AlignmentUtil();
    public static List<CigarOp> structureCigar(TranscriptFeature feature, Set<String> allowedTypes, GenomeSequenceExtractor sequenceExtractor) {
        var result = new ArrayList<CigarOp>();

        var features = feature.getFeatures()
                .stream()
                .filter(f -> allowedTypes.contains(GtfConfig.getDefault(f.getBaseData().getType())))
                .sorted(Comparator.comparingInt(f -> f.getBaseData().getStart()))
                .toList();

        for (var g : features){
            var baseData = g.getBaseData();
            var type = baseData.getType();

            var start = baseData.getStart();
            var stop = baseData.getEnd();

            var len = stop - start + 1;
            String seq = null;

            if(features.size() == 1){
                try {
                    seq = sequenceExtractor.getSequence(baseData.getContig(), start, stop);
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

                if (maxLen <= SHORT_EXON_THRESHOLD) {
                    posSimilarity = posSimilarity >= STRICT_SHORT_MIN_SIM ? posSimilarity : 0.0;
                }

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

    public static boolean isSimilar(List<CigarOp> a, List<CigarOp> b) {
        if(a.size() == 1 && b.size() == 1) {
            return handleSingleFeature(a, b) >= SEQ_HOMOLOGY_CUTOFF;
        }

        return cigarSimilarity(a, b) >=  SIMILARITY_CUTOFF;
    }
}
