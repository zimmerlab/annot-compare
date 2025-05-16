package com.github.zimmerlab.gtfcompare.model.comparison;

import java.util.ArrayList;
import java.util.List;

public class FeatureComparisonResult {
    private String featureType;
    private boolean missinInTargetTranscript;
    private boolean missingInQueryTranscript;
    private boolean areSameFeatures = true;
    private List<RegionComparisonResult> regionComparisons = new ArrayList<>();
    private List<String> messages = new ArrayList<>();

    public String getFeatureType() {
        return featureType;
    }

    public void setFeatureType(String featureType) {
        this.featureType = featureType;
    }

    public boolean isMissingInTargetTranscript() {
        return missinInTargetTranscript;
    }

    public void setMissingInTargetTranscript(boolean missinInTargetTranscript) {
        this.missinInTargetTranscript = missinInTargetTranscript;
    }

    public boolean isMissingInQueryTranscript() {
        return missingInQueryTranscript;
    }

    public void setMissingInQueryTranscript(boolean missingInQueryTranscript) {
        this.missingInQueryTranscript = missingInQueryTranscript;
    }

    public List<RegionComparisonResult> getRegionComparisons() {
        return regionComparisons;
    }

    public void addRegionComparison(RegionComparisonResult regionComparison) {
        this.regionComparisons.add(regionComparison);
    }

    public List<String> getMessages() {
        return messages;
    }

    public void addMessage(String message) {
        this.messages.add(message);
    }

    public boolean areSameFeatures() {
        return areSameFeatures;
    }

    public void setAreSameFeatures(boolean areSameFeatures) {
        this.areSameFeatures = areSameFeatures;
    }
}