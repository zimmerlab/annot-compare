package com.github.zimmerlab.gtfcompare.compare.comparators.transcriptfeatures;

import com.github.zimmerlab.gtfcompare.compare.CDSComparisonFeature;
import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.ComparisonFeature;
import com.github.zimmerlab.gtfcompare.utils.Constants;

public class StopComparator implements ComparisonFeature, CDSComparisonFeature {
    @Override
    public String getName() {
        return Constants.STOP_COMPARATOR_NAME;
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
