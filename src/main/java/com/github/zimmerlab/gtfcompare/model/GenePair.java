package com.github.zimmerlab.gtfcompare.model;

import com.github.kleinsamuel.gtfutils.feature.GeneFeature;
import com.github.zimmerlab.gtfcompare.model.comparison.GeneComparisonResult;

public class GenePair extends FeaturePair<GeneFeature> {
    private GeneComparisonResult geneComparisonResult;

    public GenePair(GeneFeature targetGene, GeneFeature queryGene) {
        super(targetGene, queryGene);
    }

    public GenePair(GeneFeature targetGene, GeneFeature queryGene, GeneComparisonResult geneComparisonResult) {
        super(targetGene, queryGene);
        this.geneComparisonResult = geneComparisonResult;
    }

    public GeneComparisonResult getGeneComparisonResult() {
        return geneComparisonResult;
    }
}
