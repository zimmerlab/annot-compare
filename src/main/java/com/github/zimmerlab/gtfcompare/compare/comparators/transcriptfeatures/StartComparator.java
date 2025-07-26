package com.github.zimmerlab.gtfcompare.compare.comparators.transcriptfeatures;

import com.github.zimmerlab.gtfcompare.compare.CDSComparisonFeature;
import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.ComparisonFeature;
import com.github.zimmerlab.gtfcompare.utils.Constants;

public class StartComparator implements ComparisonFeature, CDSComparisonFeature {
    @Override
    public String getName() {
        return Constants.START_COMPARATOR_NAME;
    }

    @Override
    public boolean compare(ComparisonContext ctx) {
        var targetStart = ctx.getTargetFeature().getBaseData().getStart();
        var queryStart  = ctx.getQueryFeature().getBaseData().getStart();

        var thr = ctx.getConfig().getThreshold(getName());

        return Math.abs(targetStart - queryStart) > thr;
    }
}
