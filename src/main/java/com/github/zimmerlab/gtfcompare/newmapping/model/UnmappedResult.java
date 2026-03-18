package com.github.zimmerlab.gtfcompare.newmapping.model;

import java.util.Map;

public record UnmappedResult(Map<String, String> unmappedQueries,
                             Map<String, String> unmappedTargets) {
}
