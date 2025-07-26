package com.github.zimmerlab.gtfcompare.compare.comparators.transcriptfeatures;

import com.github.zimmerlab.gtfcompare.compare.CDSComparisonFeature;
import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.ComparisonFeature;
import com.github.zimmerlab.gtfcompare.utils.Constants;

public class LengthComparator implements ComparisonFeature, CDSComparisonFeature {

    public LengthComparator() {}
    @Override
    public String getName() {
        return Constants.LENGTH_COMPARATOR_NAME;
    }

    @Override
    public boolean compare(ComparisonContext ctx) {
        var targetBaseData = ctx.getTargetFeature().getBaseData();
        var queryBaseData = ctx.getQueryFeature().getBaseData();

        var targetStart = targetBaseData.getStart();
        var targetStop = targetBaseData.getEnd();

        var queryStart = queryBaseData.getStart();
        var queryStop = queryBaseData.getEnd();

        var thr = ctx.getConfig().getThreshold(getName());

        return Math.abs((targetStop - targetStart) - (queryStop - queryStart)) > thr;
    }
}
