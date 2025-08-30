package com.github.zimmerlab.gtfcompare.mapping;

public class OutputLine {
    private final int orderKey;
    private final String[] columns;

    public OutputLine(int orderKey, String... columns) {
        this.orderKey = orderKey;
        this.columns  = columns;
    }

    public int getOrderKey() {
        return orderKey;
    }

    public String[] getColumns() {
        return columns;
    }
}
