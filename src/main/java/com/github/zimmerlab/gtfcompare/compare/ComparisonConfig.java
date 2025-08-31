package com.github.zimmerlab.gtfcompare.compare;

import com.github.zimmerlab.gtfcompare.model.Impact;
import com.github.zimmerlab.gtfcompare.model.config.FeatureConfig;

import java.util.Map;
import java.util.Set;

public class ComparisonConfig {
    private final Set<String> enabledFeatures;
    private final Set<String> enabledTranscriptFeatures;
    private final Set<String> allowedGeneBiotypes;
    private final Map<String, Double> thresholds;
    private final Map<String, Impact> impactLevels;

    ComparisonConfig(Set<String> enabledFeatures, Set<String> enabledTranscriptFeatures, Map<String, Double> thresholds, Set<String> allowedGeneBiotypes, Map<String, Impact> impactLevels) {
        this.enabledFeatures = enabledFeatures;
        this.enabledTranscriptFeatures = enabledTranscriptFeatures;
        this.thresholds = thresholds;
        this.allowedGeneBiotypes = allowedGeneBiotypes;
        this.impactLevels = impactLevels;
    }

    public boolean isEnabled(String feature) {
        return enabledFeatures.contains(feature) || enabledTranscriptFeatures.contains(feature);
    }

    public Double getThreshold(String feature) {
        return thresholds.getOrDefault(feature, 0.0);
    }

    public Set<String> getAllowedGeneBiotypes() {
        return allowedGeneBiotypes;
    }

    public Map<String, Impact> getImpactLevels() {
        return impactLevels;
    }
}
