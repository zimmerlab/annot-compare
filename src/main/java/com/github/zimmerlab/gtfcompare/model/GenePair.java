package com.github.zimmerlab.gtfcompare.model;

import com.github.kleinsamuel.gtfutils.feature.GeneFeature;
import com.github.zimmerlab.gtfcompare.model.comparison.GeneComparisonResult;

public class GenePair {
    private final GeneFeature targetGene;
    private final GeneFeature queryGene;
    private GeneComparisonResult geneComparisonResult;

    public GenePair(GeneFeature targetGene, GeneFeature queryGene) {
        this.targetGene = targetGene;
        this.queryGene = queryGene;
    }

    public GenePair(GeneFeature targetGene, GeneFeature queryGene, GeneComparisonResult geneComparisonResult) {
        this.targetGene = targetGene;
        this.queryGene = queryGene;
        this.geneComparisonResult = geneComparisonResult;
    }

    public GeneFeature getTargetGene() {
        return targetGene;
    }

    public GeneFeature getQueryGene() {
        return queryGene;
    }

    public GeneComparisonResult getGeneComparisonResult() {
        return geneComparisonResult;
    }
}
