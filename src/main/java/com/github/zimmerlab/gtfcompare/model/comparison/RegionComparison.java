package com.github.zimmerlab.gtfcompare.model.comparison;

import com.sun.source.tree.UsesTree;

public class RegionComparison {
    private int start1;
    private int end1;
    private int start2;
    private int end2;
    private boolean positionDifferenceFound;
    private boolean sequenceDifferenceFound;

    public RegionComparison(int start1, int end1, int start2, int end2){
        this.start1 = start1;
        this.end1 = end1;
        this.start2 = start2;
        this.end2 = end2;
    }

    public RegionComparison(int start1, int end1, int start2, int end2, boolean positionDifferenceFound) {
        this(start1, end1, start2, end2);
        this.positionDifferenceFound = positionDifferenceFound;
    }

    public RegionComparison(boolean sequenceDifferenceFound, int start1, int end1, int start2, int end2) {
        this(start1, end1, start2, end2, false);
        this.sequenceDifferenceFound = sequenceDifferenceFound;
    }

    // Getter
    public int getStart1() {
        return start1;
    }

    public int getEnd1() {
        return end1;
    }

    public int getStart2() {
        return start2;
    }

    public int getEnd2() {
        return end2;
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

    public void setEnd1(int end1) {
        this.end1 = end1;
    }

    public void setEnd2(int end2) {
        this.end2 = end2;
    }

    public void setStart1(int start1) {
        this.start1 = start1;
    }

    public void setStart2(int start2) {
        this.start2 = start2;
    }

}
