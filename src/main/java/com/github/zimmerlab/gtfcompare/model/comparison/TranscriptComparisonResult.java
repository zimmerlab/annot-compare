package com.github.zimmerlab.gtfcompare.model.comparison;

import java.util.ArrayList;
import java.util.List;

public class TranscriptComparisonResult {
    private String transcriptId;
    private boolean transcriptMissingInGene1;
    private boolean transcriptMissingInGene2;
    private boolean startDifferent;
    private boolean stopDifferent;
    private boolean lengthDifferent;
    private boolean isSequenceDifferent;
    private List<FeatureComparisonResult> featureComparisons = new ArrayList<>();
    private List<String> messages = new ArrayList<>();

    public String getTranscriptId() {
        return transcriptId;
    }
    public void setTranscriptId(String transcriptId) {
        this.transcriptId = transcriptId;
    }
    public boolean isTranscriptMissingInGene1() {
        return transcriptMissingInGene1;
    }
    public void setTranscriptMissingInGene1(boolean transcriptMissingInGene1) {
        this.transcriptMissingInGene1 = transcriptMissingInGene1;
    }
    public boolean isTranscriptMissingInGene2() {
        return transcriptMissingInGene2;
    }
    public void setTranscriptMissingInGene2(boolean transcriptMissingInGene2) {
        this.transcriptMissingInGene2 = transcriptMissingInGene2;
    }
    public boolean isStartDifferent() {
        return startDifferent;
    }
    public void setStartDifferent(boolean startDifferent) {
        this.startDifferent = startDifferent;
    }
    public boolean isStopDifferent() {
        return stopDifferent;
    }
    public void setStopDifferent(boolean stopDifferent) {
        this.stopDifferent = stopDifferent;
    }
    public List<FeatureComparisonResult> getFeatureComparisons() {
        return featureComparisons;
    }
    public void addFeatureComparison(FeatureComparisonResult featureComparison) {
        this.featureComparisons.add(featureComparison);
    }
    public List<String> getMessages() {
        return messages;
    }
    public void addMessage(String message) {
        this.messages.add(message);
    }

    public void setLengthDifferent(boolean lengthDifferent) {
        this.lengthDifferent = lengthDifferent;
    }

    public boolean isLengthDifferent() {
        return lengthDifferent;
    }

    public void setSequenceDifferent(boolean sequenceDifferent) {
        isSequenceDifferent = sequenceDifferent;
    }

    public boolean isSequenceDifferent() {
        return isSequenceDifferent;
    }
}
