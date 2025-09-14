package com.github.zimmerlab.gtfcompare.model.config;

import com.github.kleinsamuel.gtfutils.GtfConfig;
import com.github.kleinsamuel.gtfutils.GtfConstants;
import com.github.zimmerlab.gtfcompare.model.Impact;
import com.github.zimmerlab.gtfcompare.utils.Constants;

import java.util.ArrayList;
import java.util.Map;
import java.util.Set;

public class ComparisonConfig {
    private final Set<String> enabledFeatures;
    private final Set<String> enabledTranscriptFeatures;
    private final Set<String> allowedGeneBiotypes;
    private final Map<String, Double> thresholds;
    private final Map<String, Impact> impactLevels;

    ComparisonConfig(Set<String> enabledFeatures, Set<String> enabledTranscriptFeatures, Map<String, Double> thresholds, Set<String> allowedGeneBiotypes, Map<String, Impact> impactLevels) {
        this.enabledFeatures = enabledFeatures;
        this.enabledTranscriptFeatures = enabledTranscriptFeatures;
        this.thresholds = thresholds;
        this.allowedGeneBiotypes = allowedGeneBiotypes;
        this.impactLevels = impactLevels;
    }

    public boolean isEnabled(String feature) {
        return enabledFeatures.contains(feature) || enabledTranscriptFeatures.contains(feature);
    }

    public Double getThreshold(String feature) {
        return thresholds.getOrDefault(feature, 0.0);
    }

    public Set<String> getAllowedGeneBiotypes() {
        return allowedGeneBiotypes;
    }

    public Map<String, Impact> getImpactLevels() {
        return impactLevels;
    }

    public static ComparisonConfig getComparisonConfig(ConfigJSON jsonConfig) {
        var configBuilder = new ComparisonConfigBuilder();
        configBuilder.setAllowedGeneBiotypes(jsonConfig.getGeneBiotypes().getOrDefault("allowed", new ArrayList<>()));
        var configFeatures = jsonConfig.getFeatures();

        // FEATURE COMPARATORS
        enableFeatureWithThreshold(configBuilder, Constants.LENGTH_COMPARATOR_NAME, configFeatures);
        enableFeatureWithThreshold(configBuilder, Constants.START_COMPARATOR_NAME, configFeatures);
        enableFeatureWithThreshold(configBuilder, Constants.STOP_COMPARATOR_NAME, configFeatures);

        enableFeature(configBuilder, Constants.SEQUENCE_COMPARATOR_NAME, configFeatures);
        enableFeature(configBuilder, Constants.SAME_PROTEIN_COMPARATOR_NAME, configFeatures);


        // TRANSCRIPT COMPARATORS

        enableFeatureWithThreshold(configBuilder, Constants.TRANSCRIPT_LENGTH_COMPARATOR_NAME, configFeatures);
        enableFeatureWithThreshold(configBuilder, Constants.TRANSCRIPT_START_COMPARATOR_NAME, configFeatures);
        enableFeatureWithThreshold(configBuilder, Constants.TRANSCRIPT_STOP_COMPARATOR_NAME, configFeatures);
        enableFeatureWithThreshold(configBuilder, Constants.TRANSCRIPT_STRAND_COMPARATOR_NAME, configFeatures);
        enableFeatureWithThreshold(configBuilder, Constants.TRANSCRIPT_BIOTYPE_COMPARATOR_NAME, configFeatures);


        // GENE COMPARATORS

        enableFeatureWithThreshold(configBuilder, Constants.GENE_LENGTH_COMPARATOR_NAME, configFeatures);
        enableFeatureWithThreshold(configBuilder, Constants.START_COMPARATOR_NAME, configFeatures);
        enableFeatureWithThreshold(configBuilder, Constants.STOP_COMPARATOR_NAME, configFeatures);

        enableFeature(configBuilder, Constants.GENE_STRAND_COMPARATOR_NAME, configFeatures);
        enableFeature(configBuilder, Constants.GENE_CONTIG_COMPARATOR_NAME, configFeatures);

        var transcriptFeatures = jsonConfig.getTranscriptFeatures();

        for (var transcriptFeature : Constants.FEATURE_TYPES) {
            var feature = transcriptFeatures.get(transcriptFeature);
            if (feature != null && feature.isEnabled()) {
                configBuilder.enableTranscriptFeatures(transcriptFeature);
                var impactLvl = feature.getImpactLevel();
                if(impactLvl != null){
                    configBuilder.setImpactLevels(transcriptFeature, impactLvl);

                }
                /*var th = feature.getThreshold();
                if (th != null) {
                    configBuilder.setThreshold(transcriptFeature, th);
                }*/
            }
        }

        return configBuilder.build();
    }

    private static void enableFeature(ComparisonConfigBuilder configBuilder, String featureName, Map<String, FeatureConfig> featureConfig) {
        var feature = featureConfig.get(featureName);
        if (feature != null && feature.isEnabled()) {
            configBuilder.enableFeature(featureName);
        }

        if(feature.getImpactLevel() == null) return;
        configBuilder.setImpactLevels(featureName, feature.getImpactLevel() == null ? null : feature.getImpactLevel());
    }

    private static void enableFeatureWithThreshold(ComparisonConfigBuilder configBuilder, String featureName, Map<String, FeatureConfig> featureConfig) {
        var feature = featureConfig.get(featureName);
        if (feature != null && feature.isEnabled()) {
            configBuilder.enableFeature(featureName);
            var th = feature.getThreshold();
            if (th != null) {
                configBuilder.setThreshold(featureName, th);
            }
        }

        if(feature.getImpactLevel() == null) return;
        configBuilder.setImpactLevels(featureName, feature.getImpactLevel());
    }
}
