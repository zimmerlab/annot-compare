package com.github.zimmerlab.gtfcompare.utils;

import com.github.kleinsamuel.gtfutils.GtfConfig;
import org.springframework.util.StopWatch;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import static java.util.Map.entry;

public class Constants {
    public final static List<StopWatch> STOP_WATCHES = new ArrayList<>();

    // Gene Feature Comparator Names
    public final static String GENE_LENGTH_COMPARATOR_NAME = "gene_length";
    public final static String GENE_SEQUENCE_COMPARATOR_NAME = "gene_sequence";
    public final static String GENE_START_COMPARATOR_NAME = "gene_start";
    public final static String GENE_STOP_COMPARATOR_NAME = "gene_stop";
    public final static String GENE_STRAND_COMPARATOR_NAME = "gene_strand";
    public final static String GENE_CONTIG_COMPARATOR_NAME = "gene_contig";

    // Feature Comparator Names
    public final static String LENGTH_COMPARATOR_NAME = "length";
    public final static String SEQUENCE_COMPARATOR_NAME = "sequence";
    public final static String START_COMPARATOR_NAME = "start";
    public final static String STOP_COMPARATOR_NAME = "stop";
    public final static String SAME_PROTEIN_COMPARATOR_NAME = "same_protein";

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
    public final static String INTRON = GtfConfig.TYPE_INTRON_DEFAULT;
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


    public static final Map<String, Character> CODON_TABLE = Map.ofEntries(
            // Phenylalanin (F)
            entry("TTT", 'F'), entry("TTC", 'F'),
            // Leucin (L)
            entry("TTA", 'L'), entry("TTG", 'L'),
            entry("CTT", 'L'), entry("CTC", 'L'),
            entry("CTA", 'L'), entry("CTG", 'L'),
            // Isoleucin (I)
            entry("ATT", 'I'), entry("ATC", 'I'),
            entry("ATA", 'I'),
            // Methionin (M, Start)
            entry("ATG", 'M'),
            // Valin (V)
            entry("GTT", 'V'), entry("GTC", 'V'),
            entry("GTA", 'V'), entry("GTG", 'V'),
            // Serin (S)
            entry("TCT", 'S'), entry("TCC", 'S'),
            entry("TCA", 'S'), entry("TCG", 'S'),
            entry("AGT", 'S'), entry("AGC", 'S'),
            // Prolin (P)
            entry("CCT", 'P'), entry("CCC", 'P'),
            entry("CCA", 'P'), entry("CCG", 'P'),
            // Threonin (T)
            entry("ACT", 'T'), entry("ACC", 'T'),
            entry("ACA", 'T'), entry("ACG", 'T'),
            // Alanin (A)
            entry("GCT", 'A'), entry("GCC", 'A'),
            entry("GCA", 'A'), entry("GCG", 'A'),
            // Tyrosin (Y)
            entry("TAT", 'Y'), entry("TAC", 'Y'),
            // Stop
            entry("TAA", '*'), entry("TAG", '*'),
            // Histidin (H)
            entry("CAT", 'H'), entry("CAC", 'H'),
            // Glutamin (Q)
            entry("CAA", 'Q'), entry("CAG", 'Q'),
            // Asparagin (N)
            entry("AAT", 'N'), entry("AAC", 'N'),
            // Lysin (K)
            entry("AAA", 'K'), entry("AAG", 'K'),
            // Aspartat (D)
            entry("GAT", 'D'), entry("GAC", 'D'),
            // Glutamat (E)
            entry("GAA", 'E'), entry("GAG", 'E'),
            // Cystein (C)
            entry("TGT", 'C'), entry("TGC", 'C'),
            // Stop
            entry("TGA", '*'),
            // Tryptophan (W)
            entry("TGG", 'W'),
            // Arginin (R)
            entry("CGT", 'R'), entry("CGC", 'R'),
            entry("CGA", 'R'), entry("CGG", 'R'),
            entry("AGA", 'R'), entry("AGG", 'R'),
            // Glycin (G)
            entry("GGT", 'G'), entry("GGC", 'G'),
            entry("GGA", 'G'), entry("GGG", 'G')
    );
}
