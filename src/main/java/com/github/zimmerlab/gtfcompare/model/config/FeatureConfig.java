package com.github.zimmerlab.gtfcompare.model.config;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.github.zimmerlab.gtfcompare.model.Impact;

public class FeatureConfig {
    @JsonProperty("enabled")
    private boolean enabled = false;

    @JsonProperty("threshold")
    private Double threshold;

    @JsonProperty("impact_level")
    private Impact impactLevel;

    public boolean isEnabled() {
        return enabled;
    }

    public void setEnabled(boolean enabled) {
        this.enabled = enabled;
    }

    public Double getThreshold() {
        return threshold;
    }

    public void setThreshold(Double threshold) {
        this.threshold = threshold;
    }

    public Impact getImpactLevel() {
        return impactLevel;
    }

    public void setImpactLevel(Impact impactLevel) {
        this.impactLevel = impactLevel;
    }
}