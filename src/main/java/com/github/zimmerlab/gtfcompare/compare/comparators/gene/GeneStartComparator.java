package com.github.zimmerlab.gtfcompare.compare.comparators.gene;

import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.GeneComparisonFeature;
import com.github.zimmerlab.gtfcompare.compare.TranscriptComparisonFeature;
import com.github.zimmerlab.gtfcompare.utils.Constants;

public class GeneStartComparator implements GeneComparisonFeature {
    @Override
    public String getName() {
        return Constants.GENE_START_COMPARATOR_NAME;
    }

    @Override
    public boolean compare(ComparisonContext ctx) {
        var targets = ctx.getTargetTranscriptFeatures();
        var queries = ctx.getQueryTranscriptFeatures();
/*
        if(targets.isEmpty() && queries.isEmpty()) {
            return false;
        }

        if(targets.isEmpty() || queries.isEmpty())
            return true;*/

       /* var minTargetStart = targets.get().stream()
                .mapToLong(target -> target.getBaseData().getStart())
                .min()
                .orElse(Long.MAX_VALUE);

        var minQueryStart = queries.get().stream()
                .mapToLong(query -> query.getBaseData().getStart())
                .min()
                .orElse(Long.MAX_VALUE);*/

        var minTargetStart = ctx.getTargetFeature().getBaseData().getStart();
        var minQueryStart = ctx.getQueryFeature().getBaseData().getStart();
        return Math.abs(minTargetStart - minQueryStart) > ctx.getConfig().getThreshold(getName());
    }
}
