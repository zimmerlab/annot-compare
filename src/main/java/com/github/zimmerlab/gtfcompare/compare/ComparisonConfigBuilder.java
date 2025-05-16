package com.github.zimmerlab.gtfcompare.compare;

import com.github.zimmerlab.gtfcompare.model.Impact;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class ComparisonConfigBuilder {
    private final Set<String> enabledFeatures = new HashSet<>();
    private final Map<String, Double> thresholds = new HashMap<>();
    private final Map<String, Impact> impactLevels = new HashMap<>();

    public ComparisonConfigBuilder enableFeature(String name) {
        enabledFeatures.add(name);
        return this;
    }
    public ComparisonConfigBuilder setThreshold(String name, double value) {
        thresholds.put(name, value);
        return this;
    }

    public ComparisonConfigBuilder setImpactLevels(String name, Impact level) {
        impactLevels.put(name, level);
        return this;
    }
    public ComparisonConfig build() {
        return new ComparisonConfig(enabledFeatures, thresholds);
    }
}