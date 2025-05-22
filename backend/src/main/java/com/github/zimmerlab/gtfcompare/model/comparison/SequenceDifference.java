package com.github.zimmerlab.gtfcompare.model.comparison;

public class SequenceDifference {
    private int startPosition;
    private int endPosition;
    private String diffType; // z.B. "Mismatch", "Deletion", "Insertion"
    private String details; // z.B. "seq1: A, seq2: G" oder "Deletion von 12 Basen"

    // Konstruktor, Getter und Setter
    public SequenceDifference(int startPosition, int endPosition, String diffType, String details) {
        this.startPosition = startPosition;
        this.endPosition = endPosition;
        this.diffType = diffType;
        this.details = details;
    }
    public int getStartPosition() {
        return startPosition;
    }
    public int getEndPosition() {
        return endPosition;
    }
    public String getDiffType() {
        return diffType;
    }
    public String getDetails() {
        return details;
    }
}
