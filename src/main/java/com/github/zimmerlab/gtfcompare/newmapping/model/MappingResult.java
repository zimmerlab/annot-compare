package com.github.zimmerlab.gtfcompare.newmapping.model;
import java.util.List;
import java.util.Map;
import java.util.Set;

public record MappingResult(List<ResultWithOrigin> results, Map<String, String> unmappedQueries, Map<String, String> unmappedTargets) {
}
