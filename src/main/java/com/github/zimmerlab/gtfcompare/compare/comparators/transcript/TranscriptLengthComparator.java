package com.github.zimmerlab.gtfcompare.compare.comparators.transcript;

import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.TranscriptComparisonFeature;
import com.github.zimmerlab.gtfcompare.utils.Constants;

public class TranscriptLengthComparator implements TranscriptComparisonFeature {
    @Override
    public String getName() {
        return Constants.TRANSCRIPT_LENGTH_COMPARATOR_NAME;
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
        var minTargetStart = ctx.getTargetTranscriptFeaturesMin();

        var maxQueryStop = ctx.getQueryTranscriptFeaturesMax();
        var minQueryStart = ctx.getQueryTranscriptFeaturesMin();

        if (maxTargetStop == null) {
            maxTargetStop = targets.get().stream()
                    .mapToLong(target -> target.getBaseData().getEnd())
                    .max()
                    .orElse(Long.MIN_VALUE);

            ctx.setTargetTranscriptFeaturesMax(maxTargetStop);
        }


        if(maxQueryStop == null){
            maxQueryStop = queries.get().stream()
                    .mapToLong(query -> query.getBaseData()
                            .getEnd())
                    .max()
                    .orElse(Long.MIN_VALUE);

            ctx.setQueryTranscriptFeaturesMax(maxQueryStop);
        }

        if(minTargetStart == null){
            minTargetStart = targets.get().stream()
                    .mapToLong(target -> target.getBaseData().getStart())
                    .min()
                    .orElse(Long.MAX_VALUE);

            ctx.setTargetTranscriptFeaturesMin(minTargetStart);
        }


        if(minQueryStart == null){
            minQueryStart = queries.get().stream()
                    .mapToLong(query -> query.getBaseData().getStart())
                    .min()
                    .orElse(Long.MAX_VALUE);

            ctx.setQueryTranscriptFeaturesMin(minQueryStart);
        }

        return Math.abs((maxTargetStop - minTargetStart) - (maxQueryStop - minQueryStart)) > ctx.getConfig().getThreshold(getName());
    }
}
