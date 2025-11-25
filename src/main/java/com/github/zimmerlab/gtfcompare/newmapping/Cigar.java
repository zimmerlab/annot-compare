package com.github.zimmerlab.gtfcompare.newmapping;

import com.github.kleinsamuel.gtfutils.GtfConfig;
import com.github.kleinsamuel.gtfutils.feature.GeneFeature;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;

import java.util.*;
import java.util.stream.Collectors;

import static com.github.zimmerlab.gtfcompare.newmapping.MappingConstants.SHORT_EXON_THRESHOLD;
import static com.github.zimmerlab.gtfcompare.newmapping.MappingConstants.STRICT_SHORT_MIN_SIM;

public class Cigar {
    public static List<CigarOp> structureCigar(TranscriptFeature feature, Set<String> allowedTypes) {
        var result = new ArrayList<CigarOp>();

        feature.getFeatures().stream().sorted(Comparator.comparingInt(f -> f.getBaseData().getStart())).forEach(g -> {
            var baseData = g.getBaseData();
            var type = baseData.getType();
            var defaultType = GtfConfig.getDefault(type);

            if (!allowedTypes.contains(defaultType)) return;

            var start = baseData.getStart();
            var stop = baseData.getEnd();

            var len = stop - start + 1;
            result.add(new CigarOp(type, len));
        });

        return result;
    }

    public static double cigarSimilarity(List<CigarOp> a, List<CigarOp> b) {
        if (a.isEmpty() && b.isEmpty()) return 1.0;

        if (a.isEmpty() || b.isEmpty()) return 0.0;

        if (a.size() == 1 && b.size() == 1) {
            return handleSingleExon(a, b);
        }

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

    private static double handleSingleExon(List<CigarOp> a, List<CigarOp> b) {
        var opA = a.getFirst();
        var opB = b.getFirst();

        if (!opA.type().equals(opB.type())) return 0.0;

        var maxLength = Math.max(opA.length(), opB.length());

        if (maxLength == 0) return 0.0;

        var relDiff = Math.abs(opA.length() - opB.length()) / (double) maxLength;
        var baseSimilarity = 1.0 - relDiff;

        if (maxLength < SHORT_EXON_THRESHOLD) {
            return baseSimilarity > STRICT_SHORT_MIN_SIM ? baseSimilarity : 0.0;
        }

        return baseSimilarity;
    }

    public static List<List<CigarOp>> getOrBuildCigars(Map<String, List<List<CigarOp>>> cache, GeneFeature gene, Set<String> allowedTypes) {
        return cache.computeIfAbsent(gene.getGeneId(), id -> gene.getTranscripts().stream().map(t -> structureCigar(t, allowedTypes)).collect(Collectors.toList()));
    }
}
