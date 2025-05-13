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
        var targetBaseData = ctx.getTargetFeature().getBaseData();
        var queryBaseData = ctx.getQueryFeature().getBaseData();

        var targetStart = targetBaseData.getStart();
        var queryStart = queryBaseData.getStart();

        return targetStart != queryStart;
    }
}
