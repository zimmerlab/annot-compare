package com.github.zimmerlab.gtfcompare.utils;

import com.github.zimmerlab.gtfcompare.model.FidxEntry;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
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
        var lineLength = fidx.getLineLength();
        var lineLengthWithNewline = fidx.getLineLengthWithNewLine();
        var result = new StringBuilder();
        var newLineLength = lineLengthWithNewline - lineLength;
        // var realStart = fidx.getStart();
        var linesTillStart = (start) / lineLength;
        var totalLinebreaks = linesTillStart * newLineLength;

        var realStart = fidx.getStart() + start + totalLinebreaks - 1;

        var position = 0;
        raf.seek(realStart);
        var tempBuffer = new byte[Math.min(end - start + 1, lineLength)];
        raf.readFully(tempBuffer);

        for (var b : tempBuffer) {
            var value = b & 0xFF;
            if (value < 32 || value == 127) {
                break;
            }
            position++;
        }

        raf.seek(realStart);
        tempBuffer = new byte[position];

        raf.readFully(tempBuffer);

        result.append(new String(tempBuffer));
        raf.skipBytes(newLineLength);

        var bytesToRead = end - start - position - (newLineLength - 1) + 1;
        var bytesRead = 0;

        while (bytesRead < bytesToRead) {
            var currentBytesToRead = Math.min(lineLength, bytesToRead - bytesRead);
            if (tempBuffer.length != currentBytesToRead)
                tempBuffer = new byte[currentBytesToRead];

            raf.readFully(tempBuffer);
            result.append(new String(tempBuffer));
            raf.skipBytes(newLineLength);
            bytesRead += currentBytesToRead;
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
                    case 'A':
                        yield 'T';
                    case 'T':
                        yield 'A';
                    case 'C':
                        yield 'G';
                    case 'G':
                        yield 'C';
                    default:
                        throw new Exception("yur");
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
