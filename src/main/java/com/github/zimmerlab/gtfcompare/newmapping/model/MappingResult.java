package com.github.zimmerlab.gtfcompare.newmapping.model;
import java.util.List;
import java.util.Set;

public record MappingResult(List<ResultWithOrigin> results, Set<String> unmappedQueries, Set<String> unmappedTargets) {
}
