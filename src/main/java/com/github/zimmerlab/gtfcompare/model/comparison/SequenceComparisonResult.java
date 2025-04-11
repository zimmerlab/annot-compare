package com.github.zimmerlab.gtfcompare.model.comparison;

import java.util.ArrayList;
import java.util.List;

public class SequenceComparisonResult {
    private String seq1;
    private String seq2;
    private List<SequenceDifference> differences = new ArrayList<>();

    // Getter / Setter
    public String getSeq1() {
        return seq1;
    }
    public void setSeq1(String seq1) {
        this.seq1 = seq1;
    }
    public String getSeq2() {
        return seq2;
    }
    public void setSeq2(String seq2) {
        this.seq2 = seq2;
    }
    public List<SequenceDifference> getDifferences() {
        return differences;
    }
    public void addDifference(SequenceDifference difference) {
        this.differences.add(difference);
    }
}
