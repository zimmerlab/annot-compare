package com.github.zimmerlab.gtfcompare.analysis.clique;

public class AnalysisResult {
    public final String contig;
    public final int targetGenes;
    public final int queryGenes;
    public final int clusters;
    public final int exactMatches;               // unique genes with exact match
    public final int unmatchedTargetGenes;       // no cluster at all
    public final int unmatchedQueryGenes;        // no cluster at all
    public final int inClustersNoExactTarget;    // in cluster but no exact match (target)
    public final int inClustersNoExactQuery;     // in cluster but no exact match (query)

    public AnalysisResult(String contig, int targetGenes, int queryGenes, int clusters,
                          int exactMatches, int unmatchedTargetGenes, int unmatchedQueryGenes,
                          int inClustersNoExactTarget, int inClustersNoExactQuery) {
        this.contig = contig;
        this.targetGenes = targetGenes;
        this.queryGenes = queryGenes;
        this.clusters = clusters;
        this.exactMatches = exactMatches;
        this.unmatchedTargetGenes = unmatchedTargetGenes;
        this.unmatchedQueryGenes = unmatchedQueryGenes;
        this.inClustersNoExactTarget = inClustersNoExactTarget;
        this.inClustersNoExactQuery = inClustersNoExactQuery;
    }
}