package com.github.zimmerlab.gtfcompare.newmapping.model;

import java.util.List;

public record TranscriptMappingResult(List<TranscriptIdPair> matches, List<String> unmappedTargets,
                                      List<String> unmappedQueries) {}
