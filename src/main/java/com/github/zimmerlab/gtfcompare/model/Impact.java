package com.github.zimmerlab.gtfcompare.model;

public enum Impact {
    MODERATE(3), LOW(2), HIGH(4), MODIFIER(1), NONE(0);

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
