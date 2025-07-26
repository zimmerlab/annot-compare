package com.github.zimmerlab.gtfcompare.compare;


public interface ComparisonFeature {

    String getName();

    boolean compare(ComparisonContext ctx) throws Exception;
}