package com.github.zimmerlab.gtfcompare.compare.comparators.gene;

import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.GeneComparisonFeature;
import com.github.zimmerlab.gtfcompare.utils.Constants;

public class GeneStrandComparator implements GeneComparisonFeature {
    @Override
    public String getName() {
        return Constants.GENE_STRAND_COMPARATOR_NAME;
    }

    @Override
    public boolean compare(ComparisonContext ctx) {
        var queryStrand = ctx.getQueryFeature().getBaseData().isForwardStrand();
        var targetStrand = ctx.getTargetFeature().getBaseData().isForwardStrand();
        return queryStrand != targetStrand;
    }
}
