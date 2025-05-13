package com.github.zimmerlab.gtfcompare.compare.comparators;

import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.ComparisonFeature;

public class StopComparator implements ComparisonFeature {
    @Override
    public String getName() {
        return "Stop";
    }

    @Override
    public boolean compare(ComparisonContext ctx) {
        var targetBaseData = ctx.getTargetFeature().getBaseData();
        var queryBaseData = ctx.getQueryFeature().getBaseData();

        var targetStop = targetBaseData.getEnd();
        var queryStop = queryBaseData.getEnd();

        var thr = ctx.getConfig().getThreshold(getName());

        return Math.abs(targetStop - queryStop) > thr;
    }
}
