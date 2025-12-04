package com.github.zimmerlab.gtfcompare.newmapping.outpututil;

import com.github.zimmerlab.gtfcompare.newmapping.model.ResultWithOrigin;

import java.io.IOException;
import java.io.Writer;
import java.util.List;

public class MappingOutputWriter {

    public static void write(List<ResultWithOrigin> results, Writer writer) throws IOException {
        for(var res : results) {
            var genePair = res.genePair();
            writer.write(String.format("%s\t%s\t%s\n", genePair.getQuery().getGeneId(), genePair.getTarget().getGeneId(), res.origin()));
        }
    }
}
