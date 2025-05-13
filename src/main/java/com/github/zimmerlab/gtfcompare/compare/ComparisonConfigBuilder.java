package com.github.zimmerlab.gtfcompare.compare;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class ComparisonConfigBuilder {
    private final Set<String> enabledFeatures = new HashSet<>();
    private final Map<String, Double> thresholds = new HashMap<>();

    public ComparisonConfigBuilder enableFeature(String name) {
        enabledFeatures.add(name);
        return this;
    }
    public ComparisonConfigBuilder setThreshold(String name, double value) {
        thresholds.put(name, value);
        return this;
    }
    public ComparisonConfig build() {
        return new ComparisonConfig(enabledFeatures, thresholds);
    }
}