package com.github.zimmerlab.gtfcompare.model;

public enum Impact {
    NONE(0), MODIFIER(1), LOW(2), MODERATE(3), HIGH(4);

    private final int level;

    Impact(int level) {
        this.level = level;
    }

    public static Impact max(Impact a, Impact b) {
        if (a == null) return b;
        if (b == null) return a;
        return a.level >= b.level ? a : b;
    }
}
