package com.github.zimmerlab.gtfcompare.model;

import java.util.List;

public class MappingResult<M, U> {
    private final List<M> mapping;
    private final List<U> unmappedTargets;
    private final List<U> unmappedQueries;

    public MappingResult(List<M> mapping, List<U> unmappedTargets, List<U> unmappedQueries) {
        this.mapping = mapping;
        this.unmappedTargets = unmappedTargets;
        this.unmappedQueries = unmappedQueries;
    }

    public List<U> getUnmappedQueries() {
        return unmappedQueries;
    }

    public List<U> getUnmappedTargets() {
        return unmappedTargets;
    }

    public List<M> getMapping() {
        return mapping;
    }
}
