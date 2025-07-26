package com.github.zimmerlab.gtfcompare.model;

import com.github.kleinsamuel.gtfutils.feature.GtfFeature;

public class FeaturePair {
    private final GtfFeature target;
    private final GtfFeature query;

    public FeaturePair(GtfFeature target, GtfFeature query) {
        this.target = target;
        this.query  = query;
    }

    public GtfFeature getTarget() { return target; }
    public GtfFeature getQuery()  { return query;    }
}
