package com.github.zimmerlab.gtfcompare.parser;

import com.github.zimmerlab.gtfcompare.model.Mapping;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

public class GeneMappingParser {


    public static Map<Mapping, List<Mapping>> loadChainMap(String mappingTsv, float oldRelease, float newRelease) throws IOException {
        var map = new HashMap<Mapping, List<Mapping>>();
        try (BufferedReader br = new BufferedReader(new FileReader(mappingTsv))) {
            String line = br.readLine();
            while ((line = br.readLine()) != null) {
                var f = line.split("\t");
                if (f.length != 9) continue;
                var oldIdEntry = f[1];
                var newIdEntry = f[2];
                var oldVersion = Integer.parseInt(f[3]);
                var newVersion = Integer.parseInt(f[4]);
                var oldRelEntry = Float.parseFloat(f[7]);
                var newRelEntry = Float.parseFloat(f[8]);

                if (oldRelEntry > newRelease && newRelEntry > newRelease) continue;

                var newMapping = new Mapping(newIdEntry, newRelEntry, newVersion);

                if (oldIdEntry.isEmpty()) {
                    if (newRelEntry <= newRelease) {
                        var seedValue = new Mapping(newIdEntry, newRelEntry, newVersion);

                        map.computeIfAbsent(seedValue, k -> new ArrayList<>()).add(seedValue);
                    }
                    continue;
                }

                var oldKey = new Mapping(oldIdEntry, oldRelEntry, oldVersion);
                map.computeIfAbsent(oldKey, k -> new ArrayList<>());

                if (newIdEntry.isEmpty() && newRelEntry <= oldRelease) {
                    map.remove(oldKey);
                    continue;
                }

                if (newRelEntry >= oldRelEntry && newRelEntry <= newRelease) {
                    map.get(oldKey).add(newMapping);

                    if (!oldIdEntry.equals(newIdEntry) && oldRelEntry < oldRelease) {
                        map.remove(oldKey);
                    }
                }
            }
        }

        return map;
    }

    public static Map<String, List<String>> makeFinalMap(Map<Mapping, List<Mapping>> adjacency, float startRelease, float endRelease) {
        Map<String, Mapping> startNodes = new HashMap<>();
        for (Mapping key : adjacency.keySet()) {
            String gid = key.getGeneId();
            float rel = key.getRelease();
            if (rel <= startRelease) {
                Mapping prev = startNodes.get(gid);
                if (prev == null || prev.getRelease() < rel) {
                    startNodes.put(gid, key);
                }
            }
        }

        Map<String, List<String>> result = new HashMap<>();

        for (var entry : startNodes.entrySet()) {
            String oldId = entry.getKey();
            Mapping start = entry.getValue();

            Set<Mapping> visited = new HashSet<>();
            Deque<Mapping> stack = new ArrayDeque<>();
            stack.push(start);
            visited.add(start);
            while (!stack.isEmpty()) {
                Mapping cur = stack.pop();
                for (Mapping succ : adjacency.getOrDefault(cur, List.of())) {
                    if (succ.getRelease() <= endRelease && visited.add(succ)) {
                        stack.push(succ);
                    }
                }
            }

            Map<String, Mapping> bestById = new HashMap<>();
            for (Mapping m : visited) {
                float rel = m.getRelease();
                if (rel > endRelease) continue;
                String gid = m.getGeneId();
                Mapping prev = bestById.get(gid);
                if (prev == null || rel > prev.getRelease()) {
                    bestById.put(gid, m);
                }
            }

            List<String> allIds = bestById.values().stream().map(Mapping::getGeneId).distinct().collect(Collectors.toList());

            List<String> withoutOld = allIds.stream().filter(gid -> !gid.equals(oldId)).collect(Collectors.toList());

            List<String> finals = withoutOld.isEmpty() ? allIds : withoutOld;

            result.put(oldId, finals);
        }

        return result;
    }

}