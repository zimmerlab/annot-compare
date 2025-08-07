package com.github.zimmerlab.gtfcompare.model;

import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.model.comparison.TranscriptComparisonResult;

public class TranscriptPair {
    private final TranscriptFeature targetTranscript;
    private final TranscriptFeature queryTranscript;
    private TranscriptComparisonResult transcriptComparisonResult;

    public TranscriptPair(TranscriptFeature transcript1, TranscriptFeature transcript2){
        this.targetTranscript = transcript1;
        this.queryTranscript = transcript2;
    }

    public TranscriptPair(TranscriptFeature transcript1, TranscriptFeature transcript2, TranscriptComparisonResult transcriptComparisonResult){
        this.targetTranscript = transcript1;
        this.queryTranscript = transcript2;
        this.transcriptComparisonResult = transcriptComparisonResult;
    }


    public TranscriptFeature getTargetTranscript() {
        return targetTranscript;
    }

    public TranscriptFeature getQueryTranscript() {
        return queryTranscript;
    }

    public TranscriptComparisonResult getTranscriptComparisonResult() {
        return transcriptComparisonResult;
    }

}
