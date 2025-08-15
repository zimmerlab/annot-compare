package com.github.zimmerlab.gtfcompare.compare.comparators.transcript;

import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.TranscriptComparisonFeature;
import com.github.zimmerlab.gtfcompare.utils.Constants;

public class TranscriptStopComparator implements TranscriptComparisonFeature {
    @Override
    public String getName() {
        return Constants.TRANSCRIPT_STOP_COMPARATOR_NAME;
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

        var maxTargetStop = ctx.getTargetTranscriptFeaturesMax();
        var maxQueryStop = ctx.getQueryTranscriptFeaturesMax();

        if (maxTargetStop == null) {
            maxTargetStop = targets.get().stream()
                    .mapToLong(target -> target.getBaseData().getEnd())
                    .max()
                    .orElse(Long.MIN_VALUE);

            ctx.setTargetTranscriptFeaturesMax(maxTargetStop);
        }

        if (maxQueryStop == null) {
            maxQueryStop = queries.get().stream()
                    .mapToLong(query -> query.getBaseData()
                            .getEnd())
                    .max()
                    .orElse(Long.MIN_VALUE);

            ctx.setQueryTranscriptFeaturesMax(maxQueryStop);
        }

        return Math.abs(maxTargetStop - maxQueryStop) > ctx.getConfig().getThreshold(getName());
    }
}
