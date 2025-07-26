package com.github.zimmerlab.gtfcompare.model;

import java.util.Objects;

public class Mapping {
    private String geneId;
    private String oldGeneId;
    private float release;
    private int version;
    private int mappingSession;

    public Mapping(String geneId, String oldGeneId, float release, int version, int mappingSession) {
        this.geneId = geneId;
        this.oldGeneId = oldGeneId;
        this.release = release;
        this.version = version;
        this.mappingSession = mappingSession;
    }

    public float getRelease() {
        return release;
    }

    public String getGeneId() {
        return geneId;
    }

    public String getOldGeneId() {
        return oldGeneId;
    }

    public int getMappingSession() {
        return mappingSession;
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
