package com.github.zimmerlab.gtfcompare.model.comparison;

import java.util.ArrayList;
import java.util.List;

public class FeatureComparisonResult {
    private String featureType;
    private boolean missingInTranscript1;
    private boolean missingInTranscript2;
    private List<RegionComparison> regionComparisons = new ArrayList<>();
    private List<String> messages = new ArrayList<>();

    public String getFeatureType() {
        return featureType;
    }

    public void setFeatureType(String featureType) {
        this.featureType = featureType;
    }

    public boolean isMissingInTranscript1() {
        return missingInTranscript1;
    }

    public void setMissingInTranscript1(boolean missingInTranscript1) {
        this.missingInTranscript1 = missingInTranscript1;
    }

    public boolean isMissingInTranscript2() {
        return missingInTranscript2;
    }

    public void setMissingInTranscript2(boolean missingInTranscript2) {
        this.missingInTranscript2 = missingInTranscript2;
    }

    public List<RegionComparison> getRegionComparisons() {
        return regionComparisons;
    }

    public void addRegionComparison(RegionComparison regionComparison) {
        this.regionComparisons.add(regionComparison);
    }

    public List<String> getMessages() {
        return messages;
    }

    public void addMessage(String message) {
        this.messages.add(message);
    }
}