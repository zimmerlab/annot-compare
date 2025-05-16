package com.github.zimmerlab.gtfcompare.model.comparison;

import java.util.ArrayList;
import java.util.List;

public class SequenceComparisonResult {
    private String targetSeq;
    private String querySeq;
    private boolean isSameSequence = true;
    private List<SequenceDifference> differences = new ArrayList<>();

    // Getter / Setter
    public String getTargetSeq() {
        return targetSeq;
    }
    public void setTargetSeq(String targetSeq) {
        this.targetSeq = targetSeq;
    }
    public String getQuerySeq() {
        return querySeq;
    }
    public void setQuerySeq(String querySeq) {
        this.querySeq = querySeq;
    }
    public List<SequenceDifference> getDifferences() {
        return differences;
    }
    public void addDifference(SequenceDifference difference) {
        this.differences.add(difference);
    }

    public boolean isSameSequence() {
        return isSameSequence;
    }

    public void setSameSequence(boolean sameSequence) {
        isSameSequence = sameSequence;
    }
}
