package com.github.zimmerlab.gtfcompare.model.comparison;

import java.util.ArrayList;
import java.util.List;

public class ComparisonResult {
    private String geneId1;
    private String geneId2;
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

    public String getGeneId1() {
        return geneId1;
    }

    public void setGeneId1(String geneId) {
        this.geneId1 = geneId;
    }

    public String getGeneId2() {
        return geneId2;
    }

    public void setGeneId2(String geneId2) {
        this.geneId2 = geneId2;
    }

    public boolean areSameGene() {
        return areSameGene;
    }

    public void setAreSameGene(boolean areSameGene) {
        this.areSameGene = areSameGene;
    }
}
