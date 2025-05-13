package com.github.zimmerlab.gtfcompare.compare.comparators;

import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.ComparisonFeature;

public class StartComparator implements ComparisonFeature {
    @Override
    public String getName() {
        return "Start";
    }

    @Override
    public boolean compare(ComparisonContext ctx) {
        var targetStart = ctx.getTargetFeature().getBaseData().getStart();
        var queryStart  = ctx.getQueryFeature().getBaseData().getStart();

        var thr = ctx.getConfig().getThreshold(getName());

        return Math.abs(targetStart - queryStart) > thr;
    }
}
