package com.github.zimmerlab.gtfcompare.model.config;

import com.fasterxml.jackson.annotation.JsonProperty;

import java.util.HashMap;
import java.util.Map;

public class ConfigJSON {
    @JsonProperty("features")
    private Map<String, FeatureConfig> features = new HashMap<>();
    @JsonProperty("transcript_features")
    private Map<String, FeatureConfig> transcriptFeatures = new HashMap<>();

    public Map<String, FeatureConfig> getFeatures() {
        return features;
    }
    public void setFeatures(Map<String, FeatureConfig> features) {
        this.features = features;
    }

    public Map<String, FeatureConfig> getTranscriptFeatures() {
        return transcriptFeatures;
    }

    public void setTranscriptFeatures(Map<String, FeatureConfig> transcriptFeatures) {
        this.transcriptFeatures = transcriptFeatures;
    }
}