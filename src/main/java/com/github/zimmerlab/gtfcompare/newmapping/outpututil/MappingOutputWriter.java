package com.github.zimmerlab.gtfcompare.newmapping.outpututil;

import com.github.zimmerlab.gtfcompare.newmapping.model.ResultWithOrigin;

import java.io.IOException;
import java.io.Writer;
import java.util.List;
import java.util.stream.Collectors;

public class MappingOutputWriter {

    public static void write(List<ResultWithOrigin> results, String contig, Writer writer) throws IOException {
        for (var res : results) {
            var genePair = res.genePair();
            writer.write(String.format("%s\t%s\t%s\t%s\n", contig, genePair.getQuery().getGeneId(), genePair.getTarget().getGeneId(), res.origins().stream().map(Enum::name).collect(Collectors.joining(","))));
        }
    }
}
