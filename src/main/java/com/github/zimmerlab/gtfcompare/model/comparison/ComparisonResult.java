package com.github.zimmerlab.gtfcompare.model.comparison;

import java.util.ArrayList;
import java.util.List;

public class ComparisonResult {
    private String targetGeneId;
    private String queryGeneId;
    private boolean areSameGene = true;
    private GeneComparisonResult geneComparison = new GeneComparisonResult();
    private List<TranscriptComparisonResult> transcriptComparisons = new ArrayList<>();

    // Getter / Setter
    public GeneComparisonResult getGeneComparison() {
        return geneComparison;
    }
    public void setGeneComparison(GeneComparisonResult geneComparison) {
        this.geneComparison = geneComparison;
    }
    public List<TranscriptComparisonResult> getTranscriptComparisons() {
        return transcriptComparisons;
    }
    public void addTranscriptComparison(TranscriptComparisonResult transcriptComparison) {
        this.transcriptComparisons.add(transcriptComparison);
    }

    public String getTargetGeneId() {
        return targetGeneId;
    }

    public void setTargetGeneId(String geneId) {
        this.targetGeneId = geneId;
    }

    public String getQueryGeneId() {
        return queryGeneId;
    }

    public void setQueryGeneId(String queryGeneId) {
        this.queryGeneId = queryGeneId;
    }

    public boolean areSameGene() {
        return areSameGene;
    }

    public void setAreSameGene(boolean areSameGene) {
        this.areSameGene = areSameGene;
    }
}
