package com.github.zimmerlab.gtfcompare.compare;

public interface TranscriptComparisonFeature {
    String getName();

    boolean compare(ComparisonContext ctx);

}
