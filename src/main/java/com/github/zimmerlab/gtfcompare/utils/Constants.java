package com.github.zimmerlab.gtfcompare.utils;

import com.github.kleinsamuel.gtfutils.GtfConfig;
import org.springframework.util.StopWatch;

import java.util.ArrayList;
import java.util.List;

public class Constants {
    public final static List<StopWatch> STOP_WATCHES = new ArrayList<>();

    // Transcript Feature Comparator Names
    public final static String LENGTH_COMPARATOR_NAME = "length";
    public final static String SEQUENCE_COMPARATOR_NAME = "sequence";
    public final static String START_COMPARATOR_NAME = "start";
    public final static String STOP_COMPARATOR_NAME = "stop";

    // Transcript Comparator Names
    public final static String TRANSCRIPT_LENGTH_COMPARATOR_NAME = "transcript_length";
    public final static String TRANSCRIPT_SEQUENCE_COMPARATOR_NAME = "transcript_sequence";
    public final static String TRANSCRIPT_START_COMPARATOR_NAME = "transcript_start";
    public final static String TRANSCRIPT_STOP_COMPARATOR_NAME = "transcript_stop";

    // Feature Types
    public final static String GENE = GtfConfig.TYPE_GENE_DEFAULT;
    public final static String TRANSCRIPT = GtfConfig.TYPE_TRANSCRIPT_DEFAULT;
    public final static String EXON = GtfConfig.TYPE_EXON_DEFAULT;
    public final static String CDS = GtfConfig.TYPE_CDS_DEFAULT;
    public final static String UTR5 = GtfConfig.TYPE_FIVE_PRIME_UTR_DEFAULT;
    public final static String UTR3 = GtfConfig.TYPE_THREE_PRIME_UTR_DEFAULT;
    public final static String INTRON = "intron";
    public final static String START_CODON = GtfConfig.TYPE_START_CODON_DEFAULT;
    public final static String STOP_CODON = GtfConfig.TYPE_STOP_CODON_DEFAULT;
    public static final List<String> FEATURE_TYPES = List.of(
            GENE,
            TRANSCRIPT,
            EXON,
            CDS,
            UTR5,
            UTR3,
            INTRON,
            START_CODON,
            STOP_CODON
    );

}
