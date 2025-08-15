package com.github.zimmerlab.gtfcompare.compare.comparators.transcript;

import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.TranscriptComparisonFeature;
import com.github.zimmerlab.gtfcompare.utils.Constants;

public class TranscriptStartComparator implements TranscriptComparisonFeature {
    @Override
    public String getName() {
        return Constants.TRANSCRIPT_START_COMPARATOR_NAME;
    }

    @Override
    public boolean compare(ComparisonContext ctx) {
        var targets = ctx.getTargetTranscriptFeatures();
        var queries = ctx.getQueryTranscriptFeatures();

        var minQueryStart = ctx.getQueryTranscriptFeaturesMin();
        var minTargetStart = ctx.getTargetTranscriptFeaturesMin();

        if (targets.isEmpty() && queries.isEmpty()) {
            return false;
        }

        if (targets.isEmpty() || queries.isEmpty())
            return true;

        if (minTargetStart == null) {
            minTargetStart = targets.get().stream()
                    .mapToLong(target -> target.getBaseData().getStart())
                    .min()
                    .orElse(Long.MAX_VALUE);
            ctx.setTargetTranscriptFeaturesMin(minTargetStart);
        }

        if (minQueryStart == null) {
            minQueryStart = queries.get().stream()
                    .mapToLong(query -> query.getBaseData().getStart())
                    .min()
                    .orElse(Long.MAX_VALUE);
            ctx.setQueryTranscriptFeaturesMin(minQueryStart);
        }
        return Math.abs(minTargetStart - minQueryStart) > ctx.getConfig().getThreshold(getName());
    }
}
