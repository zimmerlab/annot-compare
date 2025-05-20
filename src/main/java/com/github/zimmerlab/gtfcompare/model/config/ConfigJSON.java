package com.github.zimmerlab.gtfcompare.model.config;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonSetter;
import com.github.kleinsamuel.gtfutils.GtfConfig;

import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;

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

    @JsonSetter("transcript_features")
    public void setTranscriptFeatures(Map<String, FeatureConfig> raw) {
        this.transcriptFeatures = raw.entrySet().stream()
                .collect(Collectors.toMap(
                        e -> GtfConfig.getDefault(e.getKey()),
                        Map.Entry::getValue
                ));
    }
}