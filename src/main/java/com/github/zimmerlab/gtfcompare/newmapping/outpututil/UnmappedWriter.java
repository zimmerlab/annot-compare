package com.github.zimmerlab.gtfcompare.newmapping.outpututil;
import java.io.IOException;
import java.io.Writer;
import java.util.Set;

public class UnmappedWriter {
    public static void write(Set<String> results, String contig, String origin, Writer writer) throws IOException {
        for(var res : results) {
            writer.write(String.format("%s\t%s\t%s\n",contig, res, origin));
        }
    }
}
