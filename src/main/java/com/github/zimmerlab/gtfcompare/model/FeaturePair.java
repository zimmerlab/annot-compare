package com.github.zimmerlab.gtfcompare.model;

import com.github.kleinsamuel.gtfutils.feature.GtfFeature;

public class FeaturePair<FeatureType extends GtfFeature> {
    private final FeatureType target;
    private final FeatureType query;

    public FeaturePair(FeatureType target, FeatureType query) {
        this.target = target;
        this.query = query;
    }

    public FeatureType getTarget() {
        return target;
    }

    public FeatureType getQuery() {
        return query;
    }
}
