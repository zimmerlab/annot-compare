package com.github.zimmerlab.gtfcompare.compare.comparators.transcript;

import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.GeneComparisonFeature;
import com.github.zimmerlab.gtfcompare.compare.TranscriptComparisonFeature;
import com.github.zimmerlab.gtfcompare.utils.Constants;

public class TranscriptStrandComparator implements TranscriptComparisonFeature {
    @Override
    public String getName() {
        return Constants.TRANSCRIPT_STRAND_COMPARATOR_NAME;
    }

    @Override
    public boolean compare(ComparisonContext ctx) {
        var queryStrand = ctx.getQueryFeature().getBaseData().isForwardStrand();
        var targetStrand = ctx.getTargetFeature().getBaseData().isForwardStrand();

        ctx.setQueryForwardStrand(queryStrand);
        ctx.setTargetForwardStrand(targetStrand);

        return queryStrand != targetStrand;
    }
}
