package com.github.zimmerlab.gtfcompare.compare.comparators.transcriptfeatures;

import com.github.zimmerlab.gtfcompare.compare.CDSComparisonFeature;
import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.ComparisonFeature;
import com.github.zimmerlab.gtfcompare.utils.Constants;

import java.util.Map;

public class SameProteinComparator implements CDSComparisonFeature {
    private static final Map<String, Character> CODON_TABLE = Constants.CODON_TABLE;

    @Override
    public String getName() {
        return Constants.SAME_PROTEIN_COMPARATOR_NAME;
    }

    @Override
    public boolean compare(ComparisonContext ctx) throws Exception {
        var targetBaseData = ctx.getTargetFeature().getBaseData();
        var queryBaseData = ctx.getQueryFeature().getBaseData();

        var targetStart = targetBaseData.getStart();
        var targetStop = targetBaseData.getEnd();

        var queryStart = queryBaseData.getStart();
        var queryStop = queryBaseData.getEnd();

        if ((targetStop - targetStart) != (queryStop - queryStart)) {
            return true; // Different lengths, hence different proteins
        }

        var targetExtractor = ctx.getTargetExtractor().isPresent() ? ctx.getTargetExtractor().get() : null;
        if (targetExtractor == null) {
            throw new Exception("Target extractor is not available for protein sequence extraction.");
        }
        var queryExtractor = ctx.getQueryExtractor().isPresent() ? ctx.getQueryExtractor().get() : null;
        if (queryExtractor == null) {
            throw new Exception("Query extractor is not available for protein sequence extraction.");
        }
        var querySequence = queryExtractor.getSequence(queryBaseData.getContig(), queryStart, queryStop);
        var targetSequence = targetExtractor.getSequence(targetBaseData.getContig(), targetStart, targetStop);

        if (querySequence == null || targetSequence == null) {
            throw new Exception("Failed to extract protein sequences for comparison.");
        }


        return !isSameProtein(targetSequence, querySequence);
    }

    private boolean isSameProtein(String targetProtein, String queryProtein) {
        for (int i = 0; i + 3 < targetProtein.length(); i += 3) {
            var targetCodon = targetProtein.substring(i, i + 3);
            var queryCodon = queryProtein.substring(i, i + 3);

            var targetAA = CODON_TABLE.get(targetCodon);
            var queryAA = CODON_TABLE.get(queryCodon);
            if (targetAA == null || queryAA == null) {
                throw new IllegalArgumentException("Invalid Codon: " + (targetAA == null ? targetCodon : queryCodon));
            }

            if (!targetAA.equals(queryAA)) {
                return false;
            }

            if (targetAA == '*') {
                return true;
            }
        }


        return true;
    }
}
