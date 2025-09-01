package com.github.zimmerlab.gtfcompare.model.config;

import com.github.zimmerlab.gtfcompare.model.Impact;

import java.util.*;

public class ComparisonConfigBuilder {
    private final Set<String> enabledFeatures = new HashSet<>();
    private final Set<String> enabledTranscriptFeatures = new HashSet<>();
    private final Set<String> allowedGeneBioTypes = new HashSet<>();
    private final Map<String, Double> thresholds = new HashMap<>();
    private final Map<String, Impact> impactLevels = new HashMap<>();

    public ComparisonConfigBuilder enableFeature(String name) {
        enabledFeatures.add(name);
        return this;
    }

    public ComparisonConfigBuilder enableTranscriptFeatures(String name) {
        enabledTranscriptFeatures.add(name);
        return this;
    }

    public ComparisonConfigBuilder setThreshold(String name, double value) {
        thresholds.put(name, value);
        return this;
    }

    public ComparisonConfigBuilder setAllowedGeneBiotypes(List<String> biotypes) {
        allowedGeneBioTypes.addAll(biotypes);
        return this;
    }

    public ComparisonConfigBuilder setImpactLevels(String name, Impact level) {
        impactLevels.put(name, level);
        return this;
    }
    public ComparisonConfig build() {
        return new ComparisonConfig(enabledFeatures, enabledTranscriptFeatures, thresholds, allowedGeneBioTypes, impactLevels);
    }
}