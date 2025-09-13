package com.github.zimmerlab.gtfcompare.compare.comparators.transcript;

import com.github.kleinsamuel.gtfutils.GtfConfig;
import com.github.zimmerlab.gtfcompare.compare.ComparisonContext;
import com.github.zimmerlab.gtfcompare.compare.TranscriptComparisonFeature;
import com.github.zimmerlab.gtfcompare.utils.Constants;

public class TranscriptBiotypeChanged implements TranscriptComparisonFeature {
    @Override
    public String getName() {
        return Constants.TRANSCRIPT_BIOTYPE_COMPARATOR_NAME;
    }

    @Override
    public boolean compare(ComparisonContext ctx) {
        var targetBiotype = ctx.getTargetTranscriptBiotype();
        var queryBiotype = ctx.getQueryTranscriptBiotype();

        if(targetBiotype.equals("Liftoff") || queryBiotype.equals("Liftoff")) return false;

        return !ctx.getTargetTranscriptBiotype().equals(ctx.getQueryTranscriptBiotype());
    }
}
