package com.github.zimmerlab.gtfcompare.utils;

import com.github.zimmerlab.gtfcompare.model.FidxEntry;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.Map;

public class GenomeSequenceExtractor {
    private RandomAccessFile raf;
    Map<String, FidxEntry> fidxData;

    private static final Logger logger = LoggerFactory.getLogger(GenomeSequenceExtractor.class);


    public GenomeSequenceExtractor(File fasta, Map<String, FidxEntry> fidxData) throws FileNotFoundException {
        raf = new RandomAccessFile(fasta, "r");
        this.fidxData = fidxData;
    }

    public String getSequence(String chr, int start, int end) throws IOException {
        var fidx = fidxData.get(chr);
        if (fidx == null) throw new IllegalArgumentException("Unknown contig: " + chr);
        if (start < 1 || end < start) throw new IllegalArgumentException("Invalid coordinates");

        final int lineLength = fidx.getLineLength();
        final int lineLengthWithNewline = fidx.getLineLengthWithNewLine();
        final int newLineLength = lineLengthWithNewline - lineLength;

        final int linesTillStart = (start - 1) / lineLength; // FIX 1
        final long totalLinebreaks = (long) linesTillStart * newLineLength;
        final long realStart = fidx.getStart() + (start - 1) + totalLinebreaks; // 1-basiert → 0-basiert

        var result = new StringBuilder();

        // Erstes Teilstück bis zum Zeilenende
        raf.seek(realStart);
        int firstChunkMax = Math.min(end - start + 1, lineLength - ((start - 1) % lineLength));
        byte[] buf = new byte[firstChunkMax];
        raf.readFully(buf);

        int position = 0;
        for (byte b : buf) {
            if (b == '\n' || b == '\r') break;
            position++;
        }
        if (position > 0) {
            result.append(new String(buf, 0, position, java.nio.charset.StandardCharsets.US_ASCII));
        }

        int totalLen = end - start + 1;
        int bytesToRead = totalLen - position; // FIX 2

        if (bytesToRead > 0) {
            raf.seek(raf.getFilePointer() + newLineLength);
        }

        int bytesRead = 0;
        while (bytesRead < bytesToRead) {
            int currentBytesToRead = Math.min(lineLength, bytesToRead - bytesRead);
            if (buf.length != currentBytesToRead) buf = new byte[currentBytesToRead];

            raf.readFully(buf);
            result.append(new String(buf, java.nio.charset.StandardCharsets.US_ASCII));
            bytesRead += currentBytesToRead;

            if (bytesRead < bytesToRead) {
                raf.seek(raf.getFilePointer() + newLineLength);
            }
        }

        return result.toString();
    }

    public String fetchExonSequence(String chr, int start, int end) throws IOException {
        FidxEntry fidx = fidxData.get(chr);
        int lineLen = fidx.getLineLength();
        int linePlusNewline = fidx.getLineLengthWithNewLine();
        int newlineLen = linePlusNewline - lineLen;

        // 1) Zero-based Offset in exon-sequence
        long zeroBasedStart = start - 1;       // GTF ist 1-based inclusive
        long fullLinesBefore = zeroBasedStart / lineLen;
        long byteOffsetBase = fullLinesBefore * newlineLen + zeroBasedStart;

        // 2) Absolute Dateioffset
        long filePos = fidx.getStart() + byteOffsetBase;
        raf.seek(filePos);

        // 3) Jetzt genau (end-start+1) Buchstaben (A/C/G/T/N) lesen, Newlines überspringen
        int toRead = end - start + 1;
        StringBuilder sb = new StringBuilder(toRead);
        while (toRead > 0) {
            int b = raf.read();
            if (b < 0) break;           // EOF (sollte nicht passieren)
            char c = (char) b;
            if (c == '\n' || c == '\r') {
                // skip newline bytes
                continue;
            }
            sb.append(c);
            toRead--;
        }
        return sb.toString();
    }

    public static String getReverseComplement(String gene, int start, int end) {
        var geneLength = gene.length();
        var seqLength = end - start;
        var subsequence = gene.substring((geneLength - start - seqLength), (geneLength - start));
        return getReverseComplement(subsequence);
    }

    public static String getReverseComplement(String gene) {
        int length = gene.length();
        char[] reverseComplement = new char[length];
        try {
            for (int i = 0; i < length; i++) {
                char base = gene.charAt(length - 1 - i);
                reverseComplement[i] = switch (base) {
                    case 'A' -> 'T';
                    case 'T' -> 'A';
                    case 'C' -> 'G';
                    case 'G' -> 'C';
                    default -> 'N';
                };

                //reverseComplement[i] = Constants.COMPLEMENT_MAP.get(base);
            }
            return new String(reverseComplement);
        } catch (Exception e) {
            logger.error("Error while trying to get reverse complement", e);
        }
        return null;
    }
}
