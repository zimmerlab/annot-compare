package com.github.zimmerlab.gtfcompare.parser;

import com.github.zimmerlab.gtfcompare.model.Mapping;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

public class EnsemblMappingParser {


    private final static HashSet<String> removed = new HashSet<>();

    public static Map<Mapping, List<Mapping>> adjacencyMapping(String mappingTsv, float oldRelease, float newRelease) throws IOException {
        var map = new HashMap<Mapping, List<Mapping>>();
        try (BufferedReader br = new BufferedReader(new FileReader(mappingTsv))) {
            String line = br.readLine();
            while ((line = br.readLine()) != null) {
                var f = line.split("\t");
                if (f.length != 9) continue;
                var mappingSession = Integer.parseInt(f[0]);
                var oldIdEntry = f[1];
                var newIdEntry = f[2];
                var oldVersion = Integer.parseInt(f[3]);
                var newVersion = Integer.parseInt(f[4]);
                var oldRelEntry = Float.parseFloat(f[7]);
                var newRelEntry = Float.parseFloat(f[8]);

                // if old and new entry is newer than newRelease, skip entry
                if (oldRelEntry > newRelease && newRelEntry > newRelease) continue;

                var newMapping = new Mapping(newIdEntry, oldIdEntry.isEmpty() ? null : oldIdEntry, newRelEntry, newVersion, mappingSession);

                if ((oldIdEntry.equals("ENSG00000242836") || newIdEntry.equals("ENSG00000242836"))) {
                    var a = 2;
                }

                // only add new entry if its release is older than newRelease, then go to next entry
                if (oldIdEntry.isEmpty()) {
                    if (newRelEntry <= newRelease) {
                        map.computeIfAbsent(newMapping, k -> new ArrayList<>()).add(newMapping);
                    }
                    continue;
                }

                // make new list if old entry does not exist yet
                var oldKey = new Mapping(oldIdEntry, null, oldRelEntry, oldVersion, mappingSession);
                map.computeIfAbsent(oldKey, k -> new ArrayList<>());

                // if old entry does not exist anymore within the range oldRelease-newRelease, remove it
                if (newIdEntry.isEmpty() && newRelEntry <= newRelease) {
                    //map.remove(oldKey);
                    removed.add(oldKey.getGeneId());
                    continue;
                }

                // if entry is in range than make new mapping from old entry to new entry
                if (newRelEntry >= oldRelEntry && newRelEntry <= newRelease) {

                    map.get(oldKey).add(newMapping);

                    // TODO stimmt evtl. nicht --> ENSG00000242836
                    // remove entries if they are not relevant
                    // meaning that start id changed before our oldRelease
                    if (!oldIdEntry.equals(newIdEntry) && oldRelEntry < oldRelease) {
                        map.remove(oldKey);
                    }
                }
            }
        }

        return map;
    }

    public static Map<String, List<String>> makeFinalMap(Map<Mapping, List<Mapping>> adjacency, float startRelease, float endRelease) {
        // get all start nodes --> the latest release for every gene id that is <= start release
        var startNodes = new HashMap<String, Mapping>();
        for (Mapping key : adjacency.keySet()) {
            var gid = key.getGeneId();
            var rel = key.getRelease();
            if (rel <= startRelease) {
                var prev = startNodes.get(gid);
                if (prev == null || prev.getRelease() < rel) {
                    startNodes.put(gid, key);
                }
            }
        }

        var result = new HashMap<String, List<String>>();

        // for every start node do DFS
        for (var entry : startNodes.entrySet()) {
            var oldId = entry.getKey();
            Mapping start = entry.getValue();

            var visited = new HashSet<Mapping>();
            var stack = new ArrayDeque<Mapping>();
            stack.push(start);
            visited.add(start);
            while (!stack.isEmpty()) {
                Mapping cur = stack.pop();
                for (Mapping succ : adjacency.getOrDefault(cur, List.of())) {
                    // if succ is <= endrelease, and was not seen yet --> add to stack
                    if (succ.getRelease() <= endRelease && visited.add(succ)) {
                        stack.push(succ);
                    }
                }
            }

            var bestTrans = new HashMap<String, Mapping>();
            for (var m : visited) {
                var oId = m.getOldGeneId();
                var nId = m.getGeneId();

                if (oId == null || oId.isEmpty() || oId.equals(nId)) {
                    continue;
                }

                var key = m.getMappingSession()
                        + "|" + oId
                        + "→" + nId;
                var prev = bestTrans.get(key);
                if (prev == null || m.getRelease() > prev.getRelease()) {
                    bestTrans.put(key, m);
                }
            }

            var filtered = bestTrans.values();

            var bestById = new HashMap<String, Mapping>();
            for (var m : filtered) {
                var rel = m.getRelease();
                if (rel > endRelease) continue;
                var gid = m.getGeneId();
                var prev = bestById.get(gid);
                if (prev == null || rel > prev.getRelease()) {
                    bestById.put(gid, m);
                }
            }

            var allIds = bestById.values().stream().map(Mapping::getGeneId).distinct().collect(Collectors.toList());

            // remove identity mapping if necessary
            if (removed.contains(oldId)) {
                allIds.remove(oldId);
            }

            result.put(oldId, allIds);
        }

        // genes that do not exist anymore have to be removed from all the successor lists
        for (var successors : result.values()) {
            successors.removeAll(removed);
        }

        result.entrySet().removeIf(entry -> entry.getValue().isEmpty());

        return result;
    }

}