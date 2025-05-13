package com.github.zimmerlab.gtfcompare.model.config;

import com.fasterxml.jackson.annotation.JsonProperty;

public class FeatureConfig {
    @JsonProperty("enabled")
    private boolean enabled = false;

    @JsonProperty("threshold")
    private Double threshold;

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
}