package com.github.zimmerlab.gtfcompare.newmapping;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.zimmerlab.gtfcompare.model.GenePair;
import com.github.zimmerlab.gtfcompare.newmapping.model.MappingOrigin;
import com.github.zimmerlab.gtfcompare.newmapping.model.ResultWithOrigin;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class MappingParser {
    private final static Logger logger = LogManager.getLogger(MappingParser.class);

    public static Map<String, List<String>> readMappingFile(String mappingFile) {
        var result = new HashMap<String, List<String>>();
        try (BufferedReader reader = new BufferedReader(new FileReader(mappingFile));) {
            reader.readLine(); // header
            while (reader.ready()) {
                var line = reader.readLine();
                if (line == null) break;
                var fields = line.split("\t");

                if (fields.length != MappingFileConstants.NUM_COLS) continue;

                var contig = fields[MappingFileConstants.CONTIG_IDX];

                result.computeIfAbsent(contig, k -> new ArrayList<>()).add(line);
            }

        } catch (FileNotFoundException e) {
            logger.error("Mapping file not found", e);
        } catch (IOException e) {
            logger.error(e.getMessage());
        }

        return result;
    }

    public static List<ResultWithOrigin> parseMappingFile(List<String> mappingFile, GtfFile targetGtf, GtfFile queryGtf) {
        var results = new ArrayList<ResultWithOrigin>();

        for (var mapping : mappingFile) {
            var fields = mapping.split("\t");

            if (fields.length != MappingFileConstants.NUM_COLS) continue;
            var targetGeneId = fields[MappingFileConstants.TARGET_GENE_IDX];
            var queryGeneId = fields[MappingFileConstants.QUERY_GENE_IDX];

            var targetGene = targetGtf.getGeneFeature(targetGeneId);
            var queryGene = queryGtf.getGeneFeature(queryGeneId);

            if (targetGene == null || queryGene == null) {
                logger.warn("Target gene or query gene not found");
                continue;
            }

            var genePair = new GenePair(targetGene, queryGene);

            var origins = new ArrayList<MappingOrigin>();
            int distance = -1;

            var parts = fields[MappingFileConstants.ORIGINS_IDX].split(",");
            for (var part : parts) {
                if(part.contains(":")){
                    distance = Integer.parseInt(part.split(":")[1]);
                    continue;
                }

                if(part.equals(MappingOrigin.TRANSCRIPT_ID_MAPPING.toString())) {
                    continue;
                }

                origins.add(MappingOrigin.valueOf(part));
            }

            var result = new ResultWithOrigin(genePair, origins, distance, fields[MappingFileConstants.TARGET_TRANSCRIPT_IDX], fields[MappingFileConstants.QUERY_TRANSCRIPT_IDX]);
            results.add(result);
        }

        return results;
    }
}
