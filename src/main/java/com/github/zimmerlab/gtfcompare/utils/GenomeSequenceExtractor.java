package com.github.zimmerlab.gtfcompare.utils;

import com.github.zimmerlab.gtfcompare.model.FidxEntry;
import htsjdk.samtools.util.BlockCompressedInputStream;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.Map;

public class GenomeSequenceExtractor implements Closeable {

    private static final Logger logger = LoggerFactory.getLogger(GenomeSequenceExtractor.class);

    private final RandomAccessFile raf;
    private final BlockCompressedInputStream bgzfStream;
    private final boolean isBgzf;
    private final Map<String, FidxEntry> fidxData;

    public GenomeSequenceExtractor(File fasta, Map<String, FidxEntry> fidxData) throws IOException {
        this.fidxData = fidxData;
        this.isBgzf = isBgzfFile(fasta);
        if (isBgzf) {
            bgzfStream = new BlockCompressedInputStream(fasta);
            raf = null;
            logger.info("Opened {} as BGZF-compressed FASTA", fasta.getName());
        } else {
            raf = new RandomAccessFile(fasta, "r");
            bgzfStream = null;
        }
    }

    // Checks for the BGZF magic bytes (gzip + deflate + extra-field flag).
    // Regular gzip has the same first 2 bytes but lacks the 0x04 extra-field byte,
    // so this distinguishes BGZF from both plain text and unsupported regular gzip.
    private static boolean isBgzfFile(File file) throws IOException {
        try (FileInputStream fis = new FileInputStream(file)) {
            byte[] magic = new byte[4];
            int read = fis.read(magic);
            return read >= 4
                && (magic[0] & 0xFF) == 0x1f
                && (magic[1] & 0xFF) == 0x8b
                && (magic[2] & 0xFF) == 0x08
                && (magic[3] & 0xFF) == 0x04;
        }
    }

    public String getSequence(String chr, int start, int end) throws IOException {
        return isBgzf ? getSequenceBgzf(chr, start, end) : getSequenceRaw(chr, start, end);
    }

    public String fetchExonSequence(String chr, int start, int end) throws IOException {
        return isBgzf ? fetchExonSequenceBgzf(chr, start, end) : fetchExonSequenceRaw(chr, start, end);
    }

    private String getSequenceRaw(String chr, int start, int end) throws IOException {
        var fidx = fidxData.get(chr);
        if (fidx == null) throw new IllegalArgumentException("Unknown contig: " + chr);
        if (start < 1 || end < start) throw new IllegalArgumentException("Invalid coordinates");

        final int lineLength = fidx.getLineLength();
        final int lineLengthWithNewline = fidx.getLineLengthWithNewLine();
        final int newLineLength = lineLengthWithNewline - lineLength;

        final int linesTillStart = (start - 1) / lineLength;
        final long totalLinebreaks = (long) linesTillStart * newLineLength;
        final long realStart = fidx.getStart() + (start - 1) + totalLinebreaks;

        var result = new StringBuilder();

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
        int bytesToRead = totalLen - position;

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

    private String fetchExonSequenceRaw(String chr, int start, int end) throws IOException {
        FidxEntry fidx = fidxData.get(chr);
        int lineLen = fidx.getLineLength();
        int linePlusNewline = fidx.getLineLengthWithNewLine();
        int newlineLen = linePlusNewline - lineLen;

        long zeroBasedStart = start - 1;
        long fullLinesBefore = zeroBasedStart / lineLen;
        long byteOffsetBase = fullLinesBefore * newlineLen + zeroBasedStart;

        raf.seek(fidx.getStart() + byteOffsetBase);

        int toRead = end - start + 1;
        StringBuilder sb = new StringBuilder(toRead);
        while (toRead > 0) {
            int b = raf.read();
            if (b < 0) break;
            char c = (char) b;
            if (c == '\n' || c == '\r') continue;
            sb.append(c);
            toRead--;
        }
        return sb.toString();
    }

    // For BGZF, fidx.getStart() is a virtual file pointer (upper 48 bits = compressed block
    // offset, lower 16 bits = offset within decompressed block). We seek to the contig start
    // using that virtual pointer, then skip forward through the decompressed bytes to reach
    // the requested position — arithmetic on decompressed offsets is identical to the raw case.
    private String getSequenceBgzf(String chr, int start, int end) throws IOException {
        var fidx = fidxData.get(chr);
        if (fidx == null) throw new IllegalArgumentException("Unknown contig: " + chr);
        if (start < 1 || end < start) throw new IllegalArgumentException("Invalid coordinates");

        final int lineLength = fidx.getLineLength();
        final int lineLengthWithNewline = fidx.getLineLengthWithNewLine();
        final int newLineLength = lineLengthWithNewline - lineLength;

        final int linesTillStart = (start - 1) / lineLength;
        final long decompressedOffset = (start - 1) + (long) linesTillStart * newLineLength;

        bgzfStream.seek(fidx.getStart());
        skipFully(bgzfStream, decompressedOffset);

        var result = new StringBuilder();

        int firstChunkMax = Math.min(end - start + 1, lineLength - ((start - 1) % lineLength));
        byte[] buf = new byte[firstChunkMax];
        readFully(bgzfStream, buf);

        int position = 0;
        for (byte b : buf) {
            if (b == '\n' || b == '\r') break;
            position++;
        }
        if (position > 0) {
            result.append(new String(buf, 0, position, java.nio.charset.StandardCharsets.US_ASCII));
        }

        int totalLen = end - start + 1;
        int bytesToRead = totalLen - position;

        if (bytesToRead > 0) {
            skipFully(bgzfStream, newLineLength);
        }

        int bytesRead = 0;
        while (bytesRead < bytesToRead) {
            int currentBytesToRead = Math.min(lineLength, bytesToRead - bytesRead);
            if (buf.length != currentBytesToRead) buf = new byte[currentBytesToRead];
            readFully(bgzfStream, buf);
            result.append(new String(buf, java.nio.charset.StandardCharsets.US_ASCII));
            bytesRead += currentBytesToRead;
            if (bytesRead < bytesToRead) {
                skipFully(bgzfStream, newLineLength);
            }
        }

        return result.toString();
    }

    private String fetchExonSequenceBgzf(String chr, int start, int end) throws IOException {
        FidxEntry fidx = fidxData.get(chr);
        int lineLen = fidx.getLineLength();
        int linePlusNewline = fidx.getLineLengthWithNewLine();
        int newlineLen = linePlusNewline - lineLen;

        long zeroBasedStart = start - 1;
        long fullLinesBefore = zeroBasedStart / lineLen;
        long byteOffsetBase = fullLinesBefore * newlineLen + zeroBasedStart;

        bgzfStream.seek(fidx.getStart());
        skipFully(bgzfStream, byteOffsetBase);

        int toRead = end - start + 1;
        StringBuilder sb = new StringBuilder(toRead);
        while (toRead > 0) {
            int b = bgzfStream.read();
            if (b < 0) break;
            char c = (char) b;
            if (c == '\n' || c == '\r') continue;
            sb.append(c);
            toRead--;
        }
        return sb.toString();
    }

    private static void readFully(InputStream in, byte[] buf) throws IOException {
        int offset = 0;
        int remaining = buf.length;
        while (remaining > 0) {
            int n = in.read(buf, offset, remaining);
            if (n < 0) throw new EOFException("Unexpected end of BGZF stream");
            offset += n;
            remaining -= n;
        }
    }

    private static void skipFully(InputStream in, long n) throws IOException {
        while (n > 0) {
            long skipped = in.skip(n);
            if (skipped <= 0) throw new EOFException("Unexpected end of BGZF stream during skip");
            n -= skipped;
        }
    }

    @Override
    public void close() throws IOException {
        if (raf != null) raf.close();
        if (bgzfStream != null) bgzfStream.close();
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
            }
            return new String(reverseComplement);
        } catch (Exception e) {
            logger.error("Error while trying to get reverse complement", e);
        }
        return null;
    }
}
