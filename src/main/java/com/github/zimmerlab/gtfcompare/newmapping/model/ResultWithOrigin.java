package com.github.zimmerlab.gtfcompare.newmapping.model;

import com.github.zimmerlab.gtfcompare.model.GenePair;

import java.util.List;

public record ResultWithOrigin(GenePair genePair, List<MappingOrigin> origins, String targetTranscriptId, String queryTranscriptId) {
}
