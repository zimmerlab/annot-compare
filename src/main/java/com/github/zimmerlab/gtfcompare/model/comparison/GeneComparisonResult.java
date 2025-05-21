package com.github.zimmerlab.gtfcompare.model.comparison;

import java.util.ArrayList;
import java.util.List;

public class GeneComparisonResult {
    private boolean startDifferent;
    private boolean stopDifferent;
    private boolean differentLength;
    private boolean sequenceDifferent;
    private boolean strandDifferent;
    private boolean contigDifferent;
    private boolean missingInTargetFile;
    private boolean missingInQueryFile;
    private boolean areSameGene = true;
    private SequenceComparisonResult sequenceComparison;

    // optional infos (gene name, errors, ...)
    private List<String> messages = new ArrayList<>();

    // Getter / Setter
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

    public boolean isStrandDifferent() {
        return strandDifferent;
    }

    public void setStrandDifferent(boolean strandDifferent) {
        this.strandDifferent = strandDifferent;
    }

    public boolean isDifferentLength() {
        return differentLength;
    }

    public void setLengthDifferent(boolean differentLength) {
        this.differentLength = differentLength;
    }

    public boolean isMissingInQueryFile() {
        return missingInQueryFile;
    }

    public void setMissingInQueryFile(boolean missingInQueryFile) {
        this.missingInQueryFile = missingInQueryFile;
    }

    public boolean isMissingInTargetFile() {
        return missingInTargetFile;
    }

    public void setMissingInTargetFile(boolean missingInTargetFile) {
        this.missingInTargetFile = missingInTargetFile;
    }

    public boolean isSequenceDifferent() {
        return sequenceDifferent;
    }
    public void setSequenceDifferent(boolean sequenceDifferent) {
        this.sequenceDifferent = sequenceDifferent;
    }

    public SequenceComparisonResult getSequenceComparison() {
        if (sequenceComparison == null) {
            sequenceComparison = new SequenceComparisonResult();
        }
        return sequenceComparison;
    }

    public void setSequenceComparison(SequenceComparisonResult sequenceComparison) {
        this.sequenceComparison = sequenceComparison;
    }

    public List<String> getMessages() {
        return messages;
    }

    public void addMessage(String message) {
        this.messages.add(message);
    }

    public void setAreSameGene(boolean areSameGene) {
        this.areSameGene = areSameGene;
    }

    public boolean areSameGene() {
        return areSameGene;
    }

    public boolean isContigDifferent() {
        return contigDifferent;
    }
    public void setContigDifferent(boolean contigDifferent) {
        this.contigDifferent = contigDifferent;
    }
}