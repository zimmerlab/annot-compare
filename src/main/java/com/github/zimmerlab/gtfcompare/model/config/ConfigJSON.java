package com.github.zimmerlab.gtfcompare.model.config;

import com.fasterxml.jackson.annotation.JsonProperty;

import java.util.HashMap;
import java.util.Map;

public class ConfigJSON {
    @JsonProperty("features")
    private Map<String, FeatureConfig> features = new HashMap<>();

    public Map<String, FeatureConfig> getFeatures() {
        return features;
    }
    public void setFeatures(Map<String, FeatureConfig> features) {
        this.features = features;
    }
}