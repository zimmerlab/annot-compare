package com.github.zimmerlab.gtfcompare.model.comparison;

import java.util.ArrayList;
import java.util.List;

public class TranscriptComparisonResult {
    private String queryTranscriptId;
    private String targetTranscriptId;
    private boolean areSameTranscript = true;
    private boolean transcriptMissingInTargetGene;
    private boolean transcriptMissingInQueryGene;
    private boolean startDifferent;
    private boolean stopDifferent;
    private boolean lengthDifferent;
    private boolean isSequenceDifferent;
    private List<FeatureComparisonResult> featureComparisons = new ArrayList<>();
    private List<String> messages = new ArrayList<>();

    public String getQueryTranscriptId() {
        return queryTranscriptId;
    }
    public void setQueryTranscriptId(String transcriptId) {
        this.queryTranscriptId = transcriptId;
    }

    public String getTargetTranscriptId() {
        return targetTranscriptId;
    }

    public void setTargetTranscriptId(String targetTranscriptId) {
        this.targetTranscriptId = targetTranscriptId;
    }

    public boolean isTranscriptMissingInTargetGene() {
        return transcriptMissingInTargetGene;
    }
    public void setTranscriptMissingInTargetGene(boolean transcriptMissingInTargetGene) {
        this.transcriptMissingInTargetGene = transcriptMissingInTargetGene;
    }
    public boolean isTranscriptMissingInQueryGene() {
        return transcriptMissingInQueryGene;
    }
    public void setTranscriptMissingInQueryGene(boolean transcriptMissingInQueryGene) {
        this.transcriptMissingInQueryGene = transcriptMissingInQueryGene;
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

    public boolean areSameTranscript() {
        return areSameTranscript;
    }

    public void setAreSameTranscript(boolean areSameTranscript) {
        this.areSameTranscript = areSameTranscript;
    }
}
