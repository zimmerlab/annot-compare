package com.github.zimmerlab.gtfcompare.mapping;

import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.model.MappingResult;
import com.github.zimmerlab.gtfcompare.model.TranscriptPair;
import com.github.zimmerlab.gtfcompare.model.comparison.TranscriptComparisonResult;
import com.github.zimmerlab.gtfcompare.utils.Constants;
import com.github.zimmerlab.gtfcompare.utils.GenomeSequenceExtractor;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class Minimap2Validator {
    private static final double MIN_IDENTITY = 95.0;
    private static final double MIN_COVERAGE = 80.0;
    private static boolean USE_CDS_ONLY = false;

    private static String buildCdna(TranscriptFeature tf, GenomeSequenceExtractor ex, boolean useCdsOnly) throws IOException {
        var exons = tf.getFeatures().stream().filter(f -> Constants.EXON.equals(f.getBaseData().getType())).sorted(Comparator.comparingInt(f -> f.getBaseData().getStart())).collect(Collectors.toList());

        var cdss = tf.getFeatures().stream().filter(f -> Constants.CDS.equals(f.getBaseData().getType())).sorted(Comparator.comparingInt(f -> f.getBaseData().getStart())).collect(Collectors.toList());

        var toUse = (useCdsOnly && !cdss.isEmpty()) ? cdss : exons;

        var forward = tf.getBaseData().isForwardStrand();
        if (!forward) Collections.reverse(toUse);

        var sb = new StringBuilder();
        for (GtfFeature seg : toUse) {
            var bd = seg.getBaseData();
            sb.append(ex.getSequence(bd.getContig(), bd.getStart(), bd.getEnd()));
        }
        var seq = sb.toString();
        return forward ? seq : reverseComplement(seq);
    }

    private static String reverseComplement(String seq) {
        var rc = new StringBuilder(seq.length());
        for (int i = seq.length() - 1; i >= 0; i--) {
            switch (seq.charAt(i)) {
                case 'A':
                    rc.append('T');
                    break;
                case 'T':
                    rc.append('A');
                    break;
                case 'C':
                    rc.append('G');
                    break;
                case 'G':
                    rc.append('C');
                    break;
                default:
                    rc.append('N');
                    break;
            }
        }
        return rc.toString();
    }


    private static final int MIN_MAPQ = 0;    // moderate uniqueness
    private static final int MIN_ALIGNED_BP = 0;    // absolute floor for very short queries
    private static final double MIN_ALIGNED_FRAC = 0.0;  // 50% of query length
    private static final double MAX_SOFTCLIP_FRAC = 1.0;  // up to 10%

    public static Map<String, String> parseAndFilterSam(Path samFile, Map<String, Integer> ref) throws IOException {
        var cigarOp = Pattern.compile("(\\d+)([MIDNSHP=X])");
        var validated = new HashMap<String, String>();

        try (var in = Files.newBufferedReader(samFile)) {
            String line;
            while ((line = in.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '@') continue;

                var f = line.split("\t");
                if (f.length < 11) continue;

                // 1) Keep primary only (drop secondary 0x100, supplementary 0x800)
                var flag = Integer.parseInt(f[1]);
                if ((flag & 0x900) != 0) continue;

                // 2) MAPQ filter (relaxed)
                var mapq = Integer.parseInt(f[4]);
                if (mapq < MIN_MAPQ) continue;

                String qname = f[0];
                String rname = f[2];
                if ("*".equals(rname)) continue;

                String cigar = f[5];
                String seq = f[9];
                int seqLen = seq.length();

                // NM
                int nm = 0;
                for (int i = 11; i < f.length; i++) {
                    if (f[i].startsWith("NM:i:")) {
                        nm = Integer.parseInt(f[i].substring(5));
                        break;
                    }
                }

                // Parse CIGAR
                var m = cigarOp.matcher(cigar);
                int M = 0, I = 0, D = 0, S = 0;
                while (m.find()) {
                    int len = Integer.parseInt(m.group(1));
                    switch (m.group(2).charAt(0)) {
                        case 'M':
                        case '=':
                        case 'X':
                            M += len;
                            break; // aligned on query & ref
                        case 'I':
                            I += len;
                            break;                      // insertion on query
                        case 'D':
                            D += len;
                            break;                      // deletion on query
                        case 'S':
                            S += len;
                            break;                      // soft-clipped on query
                        default:
                            break;                                 // ignore H,N,P here
                    }
                }

                var alignedQuery = M + I;

                var minAlignedBp = Math.max(MIN_ALIGNED_BP, (int) Math.round(MIN_ALIGNED_FRAC * seqLen));
                if (alignedQuery < minAlignedBp) continue;

                //  soft clip
                var softFrac = (alignedQuery + S == 0) ? 1.0 : (double) S / (alignedQuery + S);
                if (softFrac > MAX_SOFTCLIP_FRAC) continue;

                // coverage query
                var covQ = (seqLen == 0) ? 0.0 : (alignedQuery * 100.0 / seqLen);

                // coverage ref
                var refLen = ref.get(rname);
                var alignedRef = M + D; // auf der Referenz zählen Matches/Mismatches und Deletionen
                var covR = (refLen == null || refLen == 0) ? 0.0 : (alignedRef * 100.0 / refLen);

                // alternative: min(covQ, covR)
                var meanCov = (covQ + covR) / 2.0;

                // identity from nm
                var denom = M + I + D; // ausgerichtete Basen (query-seitig + Deletionen)
                var identity = (denom == 0) ? 0.0 : (100.0 * (1.0 - ((double) nm / denom)));

                if (covQ >= MIN_COVERAGE && covR >= MIN_COVERAGE && meanCov >= MIN_COVERAGE && identity >= MIN_IDENTITY) {
                    validated.put(qname, rname);
                }
            }
        }
        return validated;
    }

    private record FastaOut(Path path, Map<String, Integer> lengths) {
    }

    private static FastaOut writeFastaFile(List<TranscriptFeature> transcripts, GenomeSequenceExtractor seqExtractor, Path outFa) throws IOException {

        Map<String, Integer> refLenByName = new HashMap<>();
        try (var w = Files.newBufferedWriter(outFa)) {
            for (TranscriptFeature tf : transcripts) {
                String seq = buildCdna(tf, seqExtractor, USE_CDS_ONLY);
                String id = tf.getTranscriptId();
                w.write(">");
                w.write(id);
                w.newLine();
                w.write(seq); // one line is fine for minimap2
                w.newLine();
                refLenByName.put(id, seq.length());
            }
        }
        return new FastaOut(outFa, refLenByName);
    }


    private static void runMinimap2Streaming(Path minimap2Exe, Path refFa, List<TranscriptFeature> queries, GenomeSequenceExtractor querySeqExtractor, Path outSam, Path errLog, int threads) throws IOException, InterruptedException {
        // Use "-" to read queries from stdin.
        var cmd = List.of(minimap2Exe.toString(), "-ax", "asm10", "--secondary=no", "-t", Integer.toString(threads), refFa.toString(), "-");

        // Do NOT merge stderr into stdout. Keep SAM clean.
        var pb = new ProcessBuilder(cmd).redirectOutput(outSam.toFile()).redirectError(errLog != null ? errLog.toFile() : ProcessBuilder.Redirect.INHERIT.file());
        var proc = pb.start();

        // Stream queries to minimap2 stdin and then close to signal EOF.
        try (var os = proc.getOutputStream(); var w = new java.io.BufferedWriter(new java.io.OutputStreamWriter(os))) {
            for (TranscriptFeature tf : queries) {
                String seq = buildCdna(tf, querySeqExtractor, USE_CDS_ONLY);
                w.write(">");
                w.write(tf.getTranscriptId());
                w.newLine();
                w.write(seq);
                w.newLine();
            }
        }

        int exit = proc.waitFor();
        if (exit != 0) {
            throw new IOException("minimap2 exited with code " + exit + (errLog != null ? " (see " + errLog.toString() + ")" : ""));
        }
    }


    public static MappingResult<TranscriptPair, TranscriptFeature> validateWithMinimap2(List<TranscriptFeature> unmappedQueries, List<TranscriptFeature> unmappedTargets, GenomeSequenceExtractor targetSeqExtractor, GenomeSequenceExtractor querySeqExtractor, Path workDir, Path minimap2Exe, int threads, boolean useCDSOnly) throws IOException, InterruptedException {
        USE_CDS_ONLY = useCDSOnly;
        if (unmappedQueries.isEmpty() || unmappedTargets.isEmpty()) {
            return new MappingResult<>(List.of(), unmappedTargets, unmappedQueries);
        }

        var refFa = workDir.resolve("targets.tmp.fa");
        var samOut = workDir.resolve("minimap2.sam");
        var errLog = workDir.resolve("minimap2.stderr.log");

        // 1) Write targets to a regular FASTA (or prebuild an .mmi once and reuse).
        var ref = writeFastaFile(unmappedTargets, targetSeqExtractor, refFa).lengths;

        // 2) Run minimap2: stream queries via stdin; keep stderr separate from SAM.
        runMinimap2Streaming(minimap2Exe, refFa, unmappedQueries, querySeqExtractor, samOut, errLog, threads);

        // 3) Parse SAM and build mapping.
        var qToT = parseAndFilterSam(samOut, ref);

        var id2Target = unmappedTargets.stream().collect(Collectors.toMap(TranscriptFeature::getTranscriptId, tf -> tf));

        var result = new ArrayList<TranscriptPair>();
        for (TranscriptFeature qf : unmappedQueries) {
            var mappedId = qToT.get(qf.getTranscriptId());
            if (mappedId != null) {
                var target = id2Target.get(mappedId);
                if (target != null) {
                    result.add(new TranscriptPair(target, qf, new TranscriptComparisonResult()));
                }
            }
        }

        var mappedQueries = result.stream().map(TranscriptPair::getQueryTranscript).collect(Collectors.toSet());

        var mappedTargets = result.stream().map(TranscriptPair::getTargetTranscript).collect(Collectors.toSet());

        var stillUnmappedQueries = unmappedQueries.stream().filter(q -> !mappedQueries.contains(q)).toList();

        var stillUnmappedTargets = unmappedTargets.stream().filter(t -> !mappedTargets.contains(t)).toList();

        deleteQuietly(refFa, samOut, errLog);
        return new MappingResult<>(result, stillUnmappedTargets, stillUnmappedQueries);
    }


    private static void deleteQuietly(Path... paths) {
        for (Path p : paths) {
            if (p == null) continue;
            try {
                Files.deleteIfExists(p);
            } catch (IOException e) {
                System.err.println("WARN: failed to delete " + p + ": " + e.getMessage());
            }
        }
    }
}