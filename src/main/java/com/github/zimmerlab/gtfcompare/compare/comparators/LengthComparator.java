package com.github.zimmerlab.gtfcompare.compare.comparators;

import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.ComparisonFeature;

public class LengthComparator implements ComparisonFeature {

    public LengthComparator() {}
    @Override
    public String getName() {
        return "Length";
    }

    @Override
    public boolean compare(ComparisonContext ctx) {
        var targetBaseData = ctx.getTargetFeature().getBaseData();
        var queryBaseData = ctx.getQueryFeature().getBaseData();

        var targetStart = targetBaseData.getStart();
        var targetStop = targetBaseData.getEnd();

        var queryStart = queryBaseData.getStart();
        var queryStop = queryBaseData.getEnd();

        // Check if the lengths are equal
        return (targetStop - targetStart) != (queryStop - queryStart);
    }
}
