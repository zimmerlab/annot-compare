package com.github.zimmerlab.gtfcompare.model;

public class FidxEntry {
    private long start;
    private int lineLength;
    private int lineLengthWithNewLine;

    public FidxEntry(long start, int lineLength, int lineLengthWithNewLine) {
        this.start = start;
        this.lineLength = lineLength;
        this.lineLengthWithNewLine = lineLengthWithNewLine;
    }

    public long getStart() {
        return start;
    }

    public int getLineLength() {
        return lineLength;
    }

    public int getLineLengthWithNewLine() {
        return lineLengthWithNewLine;
    }

}
