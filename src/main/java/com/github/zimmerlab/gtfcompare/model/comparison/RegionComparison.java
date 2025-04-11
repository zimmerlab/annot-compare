package com.github.zimmerlab.gtfcompare.model.comparison;

public class RegionComparison {
    private int start1;
    private int end1;
    private int start2;
    private int end2;
    private boolean differenceFound;

    public RegionComparison(int start1, int end1, int start2, int end2, boolean differenceFound) {
        this.start1 = start1;
        this.end1 = end1;
        this.start2 = start2;
        this.end2 = end2;
        this.differenceFound = differenceFound;
    }
    // Getter
    public int getStart1() { return start1; }
    public int getEnd1() { return end1; }
    public int getStart2() { return start2; }
    public int getEnd2() { return end2; }
    public boolean isDifferenceFound() { return differenceFound; }
}
