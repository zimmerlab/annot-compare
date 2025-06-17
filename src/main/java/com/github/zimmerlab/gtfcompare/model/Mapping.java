package com.github.zimmerlab.gtfcompare.model;

import java.util.Objects;

public class Mapping {
    private String geneId;
    private float release;
    private int version;

    public Mapping(String geneId, float release, int version) {
        this.geneId = geneId;
        this.release = release;
        this.version = version;
    }

    public float getRelease() {
        return release;
    }

    public String getGeneId() {
        return geneId;
    }

    public int getVersion() {
        return version;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Mapping mapping = (Mapping) o;
        return geneId.equals(mapping.geneId);
    }

    @Override
    public int hashCode() {
        return Objects.hash(geneId);
    }
}
