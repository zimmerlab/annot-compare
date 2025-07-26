package com.github.zimmerlab.gtfcompare.compare.comparators.gene;

import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.GeneComparisonFeature;
import com.github.zimmerlab.gtfcompare.utils.Constants;

public class GeneContigComparator implements GeneComparisonFeature {
    @Override
    public String getName() {
        return Constants.GENE_CONTIG_COMPARATOR_NAME;
    }

    @Override
    public boolean compare(ComparisonContext ctx) {
        var queryContig = ctx.getQueryFeature().getBaseData().getContig();
        var targetContig = ctx.getTargetFeature().getBaseData().getContig();
        return !queryContig.equals(targetContig);
    }
}
