package com.github.zimmerlab.gtfcompare.newmapping;

import com.github.zimmerlab.gtfcompare.newmapping.model.MappingOrigin;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.EnumSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import static com.github.zimmerlab.gtfcompare.newmapping.NewMappingConstants.MAPPING_ORIGINS_COL;


public class FilterMappings {
    private final static Logger logger = LogManager.getLogger(FilterMappings.class);

    public static List<String> parseMappingFile(String mappingFile, Set<MappingOrigin> allowedMappings, boolean requireAll) {
        var mapping = new LinkedList<String>();

        try (var reader = Files.newBufferedReader(Path.of(mappingFile))) {
            reader.lines().skip(1).filter(line -> containsAllowedMapping(line, allowedMappings, requireAll)).forEach(mapping::add);

        } catch (IOException e) {
            logger.warn("Mapping file {} could not be read", mappingFile);
        }

        return mapping;
    }

    private static boolean containsAllowedMapping(String line, Set<MappingOrigin> allowedMappings, boolean requireAll) {
        var foundOrigins = EnumSet.noneOf(MappingOrigin.class);

        var splitLine = line.split("\t");
        var mappingOrigins = splitLine[MAPPING_ORIGINS_COL].split(",");

        for (var mappingOrigin : mappingOrigins) {
            var normalized = mappingOrigin.trim();

            if (normalized.contains(":")) {
                var split = normalized.split(":");
                normalized = split[0].trim();

                try {
                    if (MappingOrigin.valueOf(split[0]) == MappingOrigin.DISTANCE) {
                        if (split.length == 2 && split[1].equals("0")) {
                            normalized = MappingOrigin.OVERLAPPING.toString();
                        }
                    }
                } catch (IllegalArgumentException e) {
                    logger.warn("Unknown mapping origin: {}", mappingOrigin);
                    continue;
                }
            }

            try {
                var origin = MappingOrigin.valueOf(normalized);
                foundOrigins.add(origin);
            } catch (IllegalArgumentException e) {
                logger.warn("Unknown mapping origin: {}", mappingOrigin);
            }
        }

        if (requireAll) {
            return foundOrigins.containsAll(allowedMappings);
        } else {
            return foundOrigins.stream().anyMatch(allowedMappings::contains);
        }
    }

    public static void writeMappingFile(List<String> filteredMappings, String outputPath) {
        try (var writer = Files.newBufferedWriter(Path.of(outputPath), StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.CREATE)) {
            writer.write(NewMappingConstants.OUTPUT_HEADER);
            for (String filteredMapping : filteredMappings) {
                writer.write(filteredMapping);
                writer.newLine();
            }
        } catch (IOException e) {
            logger.error("Failed to write mapping file", e);
        }
    }

}
