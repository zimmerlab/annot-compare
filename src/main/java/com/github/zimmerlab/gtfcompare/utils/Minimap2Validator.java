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

    private static String buildCdna(TranscriptFeature tf,
                                    GenomeSequenceExtractor ex) throws IOException {
        var exons = tf.getFeatures().stream()
                .filter(f -> "exon".equals(f.getBaseData().getType()))
                .collect(Collectors.toList());
        var cdss  = tf.getFeatures().stream()
                .filter(f -> Constants.CDS.equals(f.getBaseData().getType()))
                .collect(Collectors.toList());

        if(exons.isEmpty()){
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
                case 'A': rc.append('T'); break;
                case 'T': rc.append('A'); break;
                case 'C': rc.append('G'); break;
                case 'G': rc.append('C'); break;
                default:  rc.append('N'); break;
            }
        }
        return rc.toString();
    }


    private static void createFifo(Path fifo) throws IOException, InterruptedException {
        Files.deleteIfExists(fifo);
        var p = new ProcessBuilder("mkfifo", fifo.toString()).start();
        if (p.waitFor() != 0) {
            throw new IOException("Failed to create FIFO: " + fifo);
        }
    }

    private static void writeFifoBlocking(Path fifo, List<TranscriptFeature> transcripts, GenomeSequenceExtractor seqExtractor) throws IOException {
        try (BufferedWriter w = Files.newBufferedWriter(fifo, StandardOpenOption.WRITE)) {
            for (TranscriptFeature tf : transcripts) {
                String seq = buildCdna(tf, seqExtractor);
                System.err.printf("DEBUG-FULL: %s len=%d (expected ~%d) strand=%s%n",
                        tf.getTranscriptId(), seq.length(),
                        tf.getFeatures().stream().filter(f->"exon".equals(f.getBaseData().getType()))
                                .mapToInt(f->f.getBaseData().getEnd() - f.getBaseData().getStart() +1).sum(),
                        tf.getBaseData().isForwardStrand() ? "+" : "-"
                );
                w.write(">" + tf.getTranscriptId());
                w.newLine();
                w.write(seq);
                w.newLine();
            }
        }
    }

    private static void runMinimap2OnFifos(Path minimap2Exe, Path refFifo, Path queryFifo, Path outSam, int threads) throws IOException, InterruptedException {
        var cmd = List.of(minimap2Exe.toString(), "-ax", "map-pb", "--secondary=no", "-t", Integer.toString(threads), refFifo.toString(), queryFifo.toString());
        var pb = new ProcessBuilder(cmd).redirectErrorStream(true).redirectOutput(outSam.toFile());
        var proc = pb.start();
        int exit = proc.waitFor();
        if (exit != 0) {
            throw new IOException("minimap2 exited with code " + exit);
        }
    }

    public static Map<String, String> parseAndFilterSam(Path samFile) throws IOException {
        var cigarOp = Pattern.compile("(\\d+)([MIDNSHP=X])");
        var validated = new HashMap<String, String>();

        try (var in = Files.newBufferedReader(samFile)) {
            String line;
            while ((line = in.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '@') continue;

                // split and guard: we need at least QNAME,RNAME,CIGAR,SEQLEN, plus optional tags
                String[] f = line.split("\t");
                if (f.length < 6) {
                    // malformed or unexpected line, skip it
                    continue;
                }

                String qname = f[0];
                String rname = f[2];
                if ("*".equals(rname)) continue;  // unmapped

                String cigar = f[5];
                int seqLen = f[9].length();

                // parse NM:i tag
                int nm = 0;
                for (int i = 11; i < f.length; i++) {
                    if (f[i].startsWith("NM:i:")) {
                        nm = Integer.parseInt(f[i].substring(5));
                        break;
                    }
                }

                // compute aligned length on query
                var m = cigarOp.matcher(cigar);
                int alignedLen = 0, ins = 0;
                while (m.find()) {
                    int len = Integer.parseInt(m.group(1));
                    char op  = m.group(2).charAt(0);
                    switch (op) {
                        case 'M': case '=': case 'X':
                            alignedLen += len; break;
                        case 'I':
                            ins += len;        break;
                        // D/N können für Reference‐Coverage mitgezählt werden, wenn nötig
                    }
                }
                int queryCover = alignedLen + ins;
                double coverage = queryCover * 100.0 / seqLen;
                double identity = alignedLen * 100.0 / seqLen;  // da

                if (coverage >= MIN_COVERAGE && identity >= MIN_IDENTITY) {
                    validated.put(qname, rname);
                }
            }
        }

        return validated;
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
        var refFifo = workDir.resolve("targets.fifo.fa");
        var queryFifo = workDir.resolve("queries.fifo.fa");
        var samOut = workDir.resolve("minimap2.sam");

        var executor = Executors.newFixedThreadPool(2);
        try {
            // 1) create FIFOs
            createFifo(refFifo);
            createFifo(queryFifo);

            // 2) launch writers (they block until minimap2 opens the pipes)
            executor.submit(() -> {
                try {
                    writeFifoBlocking(refFifo, unmappedTargets, targetSeqExtractor);
                } catch (IOException e) {
                    throw new UncheckedIOException(e);
                }
            });
            executor.submit(() -> {
                try {
                    writeFifoBlocking(queryFifo, unmappedQueries, querySeqExtractor);
                } catch (IOException e) {
                    throw new UncheckedIOException(e);
                }
            });

            // 3) run minimap2 reading from FIFOs
            runMinimap2OnFifos(minimap2Exe, refFifo, queryFifo, samOut, threads);

            // 4) parse SAM and filter
            var qToT = parseAndFilterSam(samOut);

            // 5) build final mapping
            var id2Target = unmappedTargets.stream().collect(Collectors.toMap(TranscriptFeature::getTranscriptId, tf -> tf));
            var result = new ArrayList<TranscriptPair>();
            for (TranscriptFeature qf : unmappedQueries) {
                var mappedId = qToT.get(qf.getTranscriptId());
                if (mappedId != null && id2Target.containsKey(mappedId)) {
                    var tPair = new TranscriptPair(id2Target.get(mappedId), qf, new TranscriptComparisonResult());
                    result.add(tPair);
                }
            }

            var mappedQueries = result.stream()
                    .map(TranscriptPair::getQueryTranscript)
                    .collect(Collectors.toSet());

            var mappedTargets = result.stream()
                    .map(TranscriptPair::getTargetTranscript)
                    .collect(Collectors.toSet());

            var stillUnmappedQueries = unmappedQueries.stream()
                    .filter(q -> !mappedQueries.contains(q))
                    .toList();

            var stillUnmappedTargets = unmappedTargets.stream()
                    .filter(t -> !mappedTargets.contains(t))
                    .toList();

            return new MappingResult<>(result, stillUnmappedTargets, stillUnmappedQueries);
        } finally {
            // shutdown writers and remove FIFOs
            executor.shutdown();
            executor.awaitTermination(5, TimeUnit.SECONDS);
            Files.deleteIfExists(refFifo);
            Files.deleteIfExists(queryFifo);
        }
    }
}