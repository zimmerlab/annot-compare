package com.github.zimmerlab.gtfcompare.compare.comparators.gene;

import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.GeneComparisonFeature;
import com.github.zimmerlab.gtfcompare.compare.TranscriptComparisonFeature;
import com.github.zimmerlab.gtfcompare.utils.Constants;

public class GeneStopComparator implements GeneComparisonFeature {
    @Override
    public String getName() {
        return Constants.GENE_STOP_COMPARATOR_NAME;
    }

    @Override
    public boolean compare(ComparisonContext ctx) {
        var targets = ctx.getTargetTranscriptFeatures();
        var queries = ctx.getQueryTranscriptFeatures();

      /*  if (targets.isEmpty() && queries.isEmpty()) {
            return false;
        }

        if (targets.isEmpty() || queries.isEmpty())
            return true;*/

    /*    var maxTargetStop = targets.get().stream()
                .mapToLong(target -> target.getBaseData().getEnd())
                .max()
                .orElse(Long.MIN_VALUE);

        var maxQueryStop = queries.get().stream()
                .mapToLong(query -> query.getBaseData()
                .getEnd())
                .max()
                .orElse(Long.MIN_VALUE);*/

        var maxTargetStop = ctx.getTargetFeature().getBaseData().getEnd();
        var maxQueryStop = ctx.getQueryFeature().getBaseData().getEnd();

        return Math.abs(maxTargetStop - maxQueryStop) > ctx.getConfig().getThreshold(getName());
    }
}
