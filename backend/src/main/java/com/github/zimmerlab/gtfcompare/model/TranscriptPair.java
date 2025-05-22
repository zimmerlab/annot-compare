package com.github.zimmerlab.gtfcompare.model;

import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.model.comparison.TranscriptComparisonResult;

public class TranscriptPair {
    private final TranscriptFeature transcript1;
    private final TranscriptFeature transcript2;
    private final TranscriptComparisonResult transcriptComparisonResult;

    public TranscriptPair(TranscriptFeature transcript1, TranscriptFeature transcript2, TranscriptComparisonResult transcriptComparisonResult){
        this.transcript1 = transcript1;
        this.transcript2 = transcript2;
        this.transcriptComparisonResult = transcriptComparisonResult;
    }

    public TranscriptFeature getTranscript1() {
        return transcript1;
    }

    public TranscriptFeature getTranscript2() {
        return transcript2;
    }

    public TranscriptComparisonResult getTranscriptComparisonResult() {
        return transcriptComparisonResult;
    }
}
