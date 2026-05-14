package com.github.zimmerlab.gtfcompare.model;

import java.util.HashMap;
import java.util.Map;

public class DisjointSet<T> {
    private final Map<T, T> parent = new HashMap<>();

    public T find(T x) {
        parent.putIfAbsent(x, x);
        if (!parent.get(x).equals(x)) {
            parent.put(x, find(parent.get(x)));
        }
        return parent.get(x);
    }

    public void union(T a, T b) {
        T ra = find(a), rb = find(b);
        if (!ra.equals(rb)) parent.put(ra, rb);
    }
}