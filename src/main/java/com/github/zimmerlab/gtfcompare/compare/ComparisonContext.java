package com.github.zimmerlab.gtfcompare.compare;

import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
import com.github.zimmerlab.gtfcompare.utils.GenomeSequenceExtractor;

import java.util.HashMap;
import java.util.List;
import java.util.Optional;

public class ComparisonContext {
    private final GtfFeature targetFeature;
    private final GtfFeature queryFeature;
    private final ComparisonConfig config;
    private final GenomeSequenceExtractor targetExtractor;
    private final GenomeSequenceExtractor queryExtractor;
    private final String targetTranscriptBiotype;
    private final String queryTranscriptBiotype;
    private final List<? extends GtfFeature> targetTranscriptFeatures;
    private final List<? extends GtfFeature> queryTranscriptFeatures;

    private Long targetTranscriptFeaturesMin;
    private Long targetTranscriptFeaturesMax;
    private Long queryTranscriptFeaturesMin;
    private Long queryTranscriptFeaturesMax;
    private Boolean lengthChanged = null;
    private Boolean isTargetForwardStrand;
    private Boolean isQueryForwardStrand;

    public ComparisonContext(GtfFeature targetFeature,
                             GtfFeature queryFeature,
                             String targetTranscriptBiotype,
                             String queryTranscriptBiotype,
                             ComparisonConfig config,
                             GenomeSequenceExtractor targetExtractor,
                             GenomeSequenceExtractor queryExtractor,
                             List<? extends GtfFeature> targetChildFeatures,
                             List<? extends GtfFeature> queryChildFeatures) {
        this.targetFeature = targetFeature;
        this.queryFeature = queryFeature;
        this.targetTranscriptBiotype = targetTranscriptBiotype;
        this.queryTranscriptBiotype = queryTranscriptBiotype;
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

    public Long getTargetTranscriptFeaturesMax() {
        return targetTranscriptFeaturesMax;
    }

    public Long getQueryTranscriptFeaturesMax() {
        return queryTranscriptFeaturesMax;
    }

    public Long getQueryTranscriptFeaturesMin() {
        return queryTranscriptFeaturesMin;
    }

    public Long getTargetTranscriptFeaturesMin() {
        return targetTranscriptFeaturesMin;
    }

    public void setQueryTranscriptFeaturesMax(Long queryTranscriptFeaturesMax) {
        this.queryTranscriptFeaturesMax = queryTranscriptFeaturesMax;
    }

    public void setQueryTranscriptFeaturesMin(Long queryTranscriptFeaturesMin) {
        this.queryTranscriptFeaturesMin = queryTranscriptFeaturesMin;
    }

    public void setTargetTranscriptFeaturesMax(Long targetTranscriptFeaturesMax) {
        this.targetTranscriptFeaturesMax = targetTranscriptFeaturesMax;
    }

    public void setTargetTranscriptFeaturesMin(Long targetTranscriptFeaturesMin) {
        this.targetTranscriptFeaturesMin = targetTranscriptFeaturesMin;
    }

    public Boolean getQueryForwardStrand() {
        return isQueryForwardStrand;
    }

    public Boolean getTargetForwardStrand() {
        return isTargetForwardStrand;
    }

    public void setQueryForwardStrand(Boolean queryForwardStrand) {
        isQueryForwardStrand = queryForwardStrand;
    }

    public void setTargetForwardStrand(Boolean targetForwardStrand) {
        isTargetForwardStrand = targetForwardStrand;
    }

    public String getQueryTranscriptBiotype() {
        return queryTranscriptBiotype;
    }

    public String getTargetTranscriptBiotype() {
        return targetTranscriptBiotype;
    }
}