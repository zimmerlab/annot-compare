package com.github.zimmerlab.gtfcompare.compare;

import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
import com.github.zimmerlab.gtfcompare.utils.GenomeSequenceExtractor;

import java.util.List;
import java.util.Optional;

public class ComparisonContext {
    private final GtfFeature targetFeature;
    private final GtfFeature queryFeature;
    private final ComparisonConfig config;
    private final GenomeSequenceExtractor targetExtractor;
    private final GenomeSequenceExtractor queryExtractor;
    private final List<? extends GtfFeature> targetTranscriptFeatures;
    private final List<? extends GtfFeature> queryTranscriptFeatures;

    private Boolean lengthChanged = null;

    public ComparisonContext(GtfFeature targetFeature,
                             GtfFeature queryFeature,
                             ComparisonConfig config,
                             GenomeSequenceExtractor targetExtractor,
                             GenomeSequenceExtractor queryExtractor,
                             List<? extends GtfFeature> targetChildFeatures,
                             List<? extends GtfFeature> queryChildFeatures) {
        this.targetFeature = targetFeature;
        this.queryFeature = queryFeature;
        this.config = config;
        this.targetExtractor = targetExtractor;
        this.queryExtractor = queryExtractor;
        this.targetTranscriptFeatures = targetChildFeatures;
        this.queryTranscriptFeatures = queryChildFeatures;
    }

    public GtfFeature getTargetFeature() {
        return targetFeature;
    }

    public GtfFeature getQueryFeature() {
        return queryFeature;
    }

    public ComparisonConfig getConfig() {
        return config;
    }

    public Optional<GenomeSequenceExtractor> getTargetExtractor() {
        return Optional.ofNullable(targetExtractor);
    }

    public Optional<GenomeSequenceExtractor> getQueryExtractor() {
        return Optional.ofNullable(queryExtractor);
    }

    public Optional<List<? extends GtfFeature>> getTargetTranscriptFeatures() {
        return Optional.ofNullable(targetTranscriptFeatures);
    }

    public Optional<List<? extends GtfFeature>> getQueryTranscriptFeatures() {
        return Optional.ofNullable(queryTranscriptFeatures);
    }

    public void setLengthChanged(Boolean lengthChanged) {
        this.lengthChanged = lengthChanged;
    }

    public Boolean hasLengthChanged() {
        return lengthChanged;
    }
}