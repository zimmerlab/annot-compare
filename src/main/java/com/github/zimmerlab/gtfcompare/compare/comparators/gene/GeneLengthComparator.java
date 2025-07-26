package com.github.zimmerlab.gtfcompare.compare.comparators.gene;

import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.GeneComparisonFeature;
import com.github.zimmerlab.gtfcompare.utils.Constants;

public class GeneLengthComparator implements GeneComparisonFeature {
    @Override
    public String getName() {
        return Constants.GENE_LENGTH_COMPARATOR_NAME;
    }

    @Override
    public boolean compare(ComparisonContext ctx) {
        var targets = ctx.getTargetTranscriptFeatures();
        var queries = ctx.getQueryTranscriptFeatures();

        if (targets.isEmpty() && queries.isEmpty()) {
            return false;
        }

        if (targets.isEmpty() || queries.isEmpty())
            return true;

        var maxTargetStop = targets.get().stream()
                .mapToLong(target -> target.getBaseData().getEnd())
                .max()
                .orElse(Long.MIN_VALUE);

        var maxQueryStop = queries.get().stream()
                .mapToLong(query -> query.getBaseData()
                .getEnd())
                .max()
                .orElse(Long.MIN_VALUE);

        var minTargetStart = targets.get().stream()
                .mapToLong(target -> target.getBaseData().getStart())
                .min()
                .orElse(Long.MAX_VALUE);

        var minQueryStart = queries.get().stream()
                .mapToLong(query -> query.getBaseData().getStart())
                .min()
                .orElse(Long.MAX_VALUE);

        return Math.abs((maxTargetStop - minTargetStart) - (maxQueryStop - minQueryStart)) > ctx.getConfig().getThreshold(getName());
    }
}
