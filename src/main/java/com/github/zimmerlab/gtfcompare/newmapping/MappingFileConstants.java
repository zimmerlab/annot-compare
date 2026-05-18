package com.github.zimmerlab.gtfcompare.newmapping;

public class MappingFileConstants {
    public static final byte NUM_COLS = 6;


    public static final byte CONTIG_IDX = 0;
    public static final byte QUERY_GENE_IDX = 1;
    public static final byte TARGET_GENE_IDX = 2;
    public static final byte QUERY_TRANSCRIPT_IDX = 3;
    public static final byte TARGET_TRANSCRIPT_IDX = 4;
    public static final byte ORIGINS_IDX = 5;

    public static final String OUTPUT_HEADER = "contig\tqueryId\ttargetId\tqueryTranscriptId\ttargetTranscriptId\tmapping_origins\n";
}
