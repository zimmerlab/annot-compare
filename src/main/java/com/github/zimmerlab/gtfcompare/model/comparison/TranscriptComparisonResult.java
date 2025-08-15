package com.github.zimmerlab.gtfcompare.model.comparison;

import java.util.ArrayList;
import java.util.List;

public class TranscriptComparisonResult {
    private String queryTranscriptId;
    private String targetTranscriptId;
    private String targetBiotype;
    private String queryBiotype;
    private String contig;
    private boolean biotypeDifferent;
    private boolean areSameTranscript = true;
    private boolean transcriptMissingInTargetGene;
    private boolean transcriptMissingInQueryGene;
    private boolean startDifferent;
    private boolean stopDifferent;
    private boolean lengthDifferent;
    private boolean isSequenceDifferent;
    private boolean isStrandDifferent;
    private Boolean isTargetForwardStrand;
    private Boolean isQueryForwardStrand;
    private List<FeatureComparisonResult> featureComparisons = new ArrayList<>();
    private List<String> messages = new ArrayList<>();
    private Long queryStart;
    private Long queryStop;
    private Long targetStart;
    private Long targetStop;

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

    public void setTargetBiotype(String targetBiotype) {
        this.targetBiotype = targetBiotype;
    }

    public void setQueryBiotype(String queryBiotype) {
        this.queryBiotype = queryBiotype;
    }

    public String getQueryBiotype() {
        return queryBiotype;
    }

    public String getTargetBiotype() {
        return targetBiotype;
    }

    public Long getQueryStart() {
        return queryStart;
    }

    public void setQueryStart(Long queryStart) {
        this.queryStart = queryStart;
    }

    public Long getQueryStop() {
        return queryStop;
    }

    public void setQueryStop(Long queryStop) {
        this.queryStop = queryStop;
    }

    public Long getTargetStart() {
        return targetStart;
    }

    public void setTargetStart(Long targetStart) {
        this.targetStart = targetStart;
    }

    public Long getTargetStop() {
        return targetStop;
    }

    public void setTargetStop(Long targetStop) {
        this.targetStop = targetStop;
    }

    public String getContig() {
        return contig;
    }

    public void setContig(String contig) {
        this.contig = contig;
    }

    public Boolean isTargetForwardStrand() {
        return isTargetForwardStrand;
    }

    public void setTargetForwardStrand(Boolean isForwardStrand) {
        this.isTargetForwardStrand = isForwardStrand;
    }

    public Boolean isQueryForwardStrand() {
        return isQueryForwardStrand;
    }

    public void setQueryForwardStrand(Boolean isForwardStrand) {
        this.isQueryForwardStrand = isForwardStrand;
    }

    public void setStrandDifferent(boolean strandDifferent) {
        isStrandDifferent = strandDifferent;
    }

    public boolean isStrandDifferent() {
        return isStrandDifferent;
    }

    public void setBiotypeDifferent(boolean biotypeDifferent) {
        this.biotypeDifferent = biotypeDifferent;
    }

    public boolean isBiotypeDifferent() {
        return biotypeDifferent;
    }
}
