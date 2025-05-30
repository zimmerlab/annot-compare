package com.github.zimmerlab.gtfcompare.model.comparison;

public class RegionComparisonResult {
    private int targetStart;
    private int targetEnd;
    private int queryStart;
    private int queryEnd;

    private boolean startDifferent;
    private boolean endDifferent;
    private boolean positionDifferenceFound;
    private boolean lengthDifferenceFound;
    private boolean sequenceDifferenceFound;
    private boolean isMissingInTargetFile;
    private boolean isMissingInQueryFile;
    private boolean isProteinDifferent;

    public RegionComparisonResult(int targetStart, int targetEnd, int queryStart, int queryEnd){
        this.targetStart = targetStart;
        this.targetEnd = targetEnd;
        this.queryStart = queryStart;
        this.queryEnd = queryEnd;
    }

    public RegionComparisonResult(int targetStart, int targetEnd, int queryStart, int queryEnd, boolean positionDifferenceFound) {
        this(targetStart, targetEnd, queryStart, queryEnd);
        this.positionDifferenceFound = positionDifferenceFound;
    }

    // Getter
    public int getTargetStart() {
        return targetStart;
    }

    public int getTargetEnd() {
        return targetEnd;
    }

    public int getQueryStart() {
        return queryStart;
    }

    public int getQueryEnd() {
        return queryEnd;
    }

    public boolean isPositionDifferenceFound() {
        return positionDifferenceFound;
    }

    public boolean isSequenceDifferenceFound() {
        return sequenceDifferenceFound;
    }

    public void setSequenceDifferenceFound(boolean sequenceDifferenceFound) {
        this.sequenceDifferenceFound = sequenceDifferenceFound;
    }

    public void setPositionDifferenceFound(boolean positionDifferenceFound) {
        this.positionDifferenceFound = positionDifferenceFound;
    }

    public void setTargetEnd(int targetEnd) {
        this.targetEnd = targetEnd;
    }

    public void setQueryEnd(int queryEnd) {
        this.queryEnd = queryEnd;
    }

    public void setTargetStart(int targetStart) {
        this.targetStart = targetStart;
    }

    public void setQueryStart(int queryStart) {
        this.queryStart = queryStart;
    }

    public boolean isLengthDifferenceFound() {
        return lengthDifferenceFound;
    }

    public void setLengthDifferent(boolean lengthDifferenceFound) {
        this.lengthDifferenceFound = lengthDifferenceFound;
    }

    public boolean isMissingInTargetFile() {
        return isMissingInTargetFile;
    }

    public boolean isMissingInQueryFile() {
        return isMissingInQueryFile;
    }

    public void setMissingInTargetFile(boolean missingInTargetFile) {
        isMissingInTargetFile = missingInTargetFile;
    }

    public void setMissingInQueryFile(boolean missingInQueryFile) {
        isMissingInQueryFile = missingInQueryFile;
    }

    public boolean isEndDifferent() {
        return endDifferent;
    }

    public boolean isStartDifferent() {
        return startDifferent;
    }

    public void setEndDifferent(boolean endDifferent) {
        this.endDifferent = endDifferent;
    }

    public void setStartDifferent(boolean startDifferent) {
        this.startDifferent = startDifferent;
    }

    public boolean isProteinDifferent() {
        return isProteinDifferent;
    }

    public void setProteinDifferent(boolean proteinDifferent) {
        isProteinDifferent = proteinDifferent;
    }
}
