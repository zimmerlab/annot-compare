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
            var origins = res.origins().stream().map(Enum::name).collect(Collectors.joining(","));
            var distance = res.geneDistance();


            writer.write(String.format("%s\t%s\t%s\t%s\t%s\t%s,DISTANCE:%s\n", contig, genePair.getQuery().getGeneId(), genePair.getTarget().getGeneId(), res.queryTranscriptId(), res.targetTranscriptId(), origins, distance));
        }
    }
}
