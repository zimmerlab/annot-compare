package com.github.zimmerlab.gtfcompare.newmapping.outpututil;
import java.io.IOException;
import java.io.Writer;
import java.util.Map;
import java.util.Set;

public class UnmappedWriter {
    public static void write(Map<String, String> results, String contig, String origin, Writer writer) throws IOException {
        for(var res : results.entrySet()) {
            writer.write(String.format("%s\t%s\t%s\t%s\n",contig, res.getKey(),res.getValue() != null ? res.getValue() : "", origin));
        }
    }
}
