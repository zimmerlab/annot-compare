package com.github.zimmerlab.gtfcompare.compare;

public interface GeneComparisonFeature {
    String getName();

    boolean compare(ComparisonContext ctx);

}
