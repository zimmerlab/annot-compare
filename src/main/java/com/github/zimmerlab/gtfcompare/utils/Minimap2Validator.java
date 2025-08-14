package com.github.zimmerlab.gtfcompare.utils;

import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.model.MappingResult;
import com.github.zimmerlab.gtfcompare.model.TranscriptPair;
import com.github.zimmerlab.gtfcompare.model.comparison.TranscriptComparisonResult;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class Minimap2Validator {
    private static final double MIN_IDENTITY = 95.0;      // percent
    private static final double MIN_COVERAGE = 80.0;      // percent

    private static String buildCdna(TranscriptFeature tf, GenomeSequenceExtractor ex) throws IOException {
        var exons = tf.getFeatures().stream().filter(f -> "exon".equals(f.getBaseData().getType())).collect(Collectors.toList());
        var cdss = tf.getFeatures().stream().filter(f -> Constants.CDS.equals(f.getBaseData().getType())).collect(Collectors.toList());

        if (exons.isEmpty()) {
            var a = 2;
        }
        var toUse = cdss.isEmpty() ? exons : cdss;

        toUse.sort(Comparator.comparingInt(f -> f.getBaseData().getStart()));
        var forward = tf.getBaseData().isForwardStrand();
        if (!forward) {
            Collections.reverse(toUse);
        }

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


    public static Map<String, String> parseAndFilterSam(Path samFile) throws IOException {
        var cigarOp = Pattern.compile("(\\d+)([MIDNSHP=X])");
        var validated = new HashMap<String, String>();
        try (var in = Files.newBufferedReader(samFile)) {
            String line;
            while ((line = in.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '@') continue;
                String[] f = line.split("\t");
                if (f.length < 11) continue; // need at least up to optional tags

                String qname = f[0];
                String rname = f[2];
                if ("*".equals(rname)) continue;

                String cigar = f[5];
                String seq = f[9];
                int seqLen = seq.length();

                int nm = 0;
                for (int i = 11; i < f.length; i++) {
                    if (f[i].startsWith("NM:i:")) {
                        nm = Integer.parseInt(f[i].substring(5));
                        break;
                    }
                }

                var m = cigarOp.matcher(cigar);
                int M = 0, I = 0, D = 0;
                while (m.find()) {
                    int len = Integer.parseInt(m.group(1));
                    switch (m.group(2).charAt(0)) {
                        case 'M':
                        case '=':
                        case 'X':
                            M += len;
                            break;
                        case 'I':
                            I += len;
                            break;
                        case 'D':
                            D += len;
                            break;
                        default:
                            break; // ignore S,H,N,P for identity here
                    }
                }

                double queryAligned = M + I;
                double coverage = (seqLen == 0) ? 0.0 : (queryAligned * 100.0 / seqLen);

                double denom = M + I + D;
                double identity = (denom == 0) ? 0.0 : (100.0 * (1.0 - ((double) nm / denom)));

                if (coverage >= MIN_COVERAGE && identity >= MIN_IDENTITY) {
                    validated.put(qname, rname);
                }
            }
        }
        return validated;
    }

    private static Path writeFastaFile(List<TranscriptFeature> transcripts, GenomeSequenceExtractor seqExtractor, Path outFa) throws IOException {
        try (var w = Files.newBufferedWriter(outFa)) {
            for (TranscriptFeature tf : transcripts) {
                String seq = buildCdna(tf, seqExtractor);
                w.write(">");
                w.write(tf.getTranscriptId());
                w.newLine();
                w.write(seq); // one line is fine for minimap2
                w.newLine();
            }
        }
        return outFa;
    }

    private static void runMinimap2Streaming(Path minimap2Exe, Path refFaOrMmi, List<TranscriptFeature> queries, GenomeSequenceExtractor querySeqExtractor, Path outSam, Path errLog, int threads) throws IOException, InterruptedException {
        // Use "-" to read queries from stdin.
        var cmd = List.of(minimap2Exe.toString(), "-ax", "map-pb", "--secondary=no", "-t", Integer.toString(threads), refFaOrMmi.toString(), "-" );

        // Do NOT merge stderr into stdout. Keep SAM clean.
        var pb = new ProcessBuilder(cmd).redirectOutput(outSam.toFile()).redirectError(errLog != null ? errLog.toFile() : ProcessBuilder.Redirect.INHERIT.file());
        var proc = pb.start();

        // Stream queries to minimap2 stdin and then close to signal EOF.
        try (var os = proc.getOutputStream(); var w = new java.io.BufferedWriter(new java.io.OutputStreamWriter(os))) {
            for (TranscriptFeature tf : queries) {
                String seq = buildCdna(tf, querySeqExtractor);
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


    /**
     * High-level driver:
     * - Creates two FIFOs
     * - Spawns writers via ExecutorService
     * - Runs minimap2
     * - Cleans up FIFOs even on error
     * - Parses SAM and returns mapping
     */
    public static MappingResult<TranscriptPair, TranscriptFeature> validateWithMinimap2(List<TranscriptFeature> unmappedQueries, List<TranscriptFeature> unmappedTargets, GenomeSequenceExtractor targetSeqExtractor, GenomeSequenceExtractor querySeqExtractor, Path workDir, Path minimap2Exe, int threads) throws IOException, InterruptedException {

        // Short-circuit when there is nothing to do.
        if (unmappedQueries.isEmpty() || unmappedTargets.isEmpty()) {
            return new MappingResult<>(List.of(), unmappedTargets, unmappedQueries);
        }

        Path refFa = workDir.resolve("targets.tmp.fa");
        Path samOut = workDir.resolve("minimap2.sam");
        Path errLog = workDir.resolve("minimap2.stderr.log");

        // 1) Write targets to a regular FASTA (or prebuild an .mmi once and reuse).
        writeFastaFile(unmappedTargets, targetSeqExtractor, refFa);

        // 2) Run minimap2: stream queries via stdin; keep stderr separate from SAM.
        runMinimap2Streaming(minimap2Exe, refFa, unmappedQueries, querySeqExtractor, samOut, errLog, threads);

        // 3) Parse SAM and build mapping.
        var qToT = parseAndFilterSam(samOut);

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