package com.github.zimmerlab.gtfcompare.model;

import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.model.comparison.TranscriptComparisonResult;

public class TranscriptPair extends FeaturePair<TranscriptFeature> {
    private TranscriptComparisonResult transcriptComparisonResult;

    public TranscriptPair(TranscriptFeature transcript1, TranscriptFeature transcript2){
        super(transcript1, transcript2);
    }

    public TranscriptPair(TranscriptFeature transcript1, TranscriptFeature transcript2, TranscriptComparisonResult transcriptComparisonResult){
        super(transcript1, transcript2);
        this.transcriptComparisonResult = transcriptComparisonResult;
    }
    
    public TranscriptComparisonResult getTranscriptComparisonResult() {
        return transcriptComparisonResult;
    }

}
