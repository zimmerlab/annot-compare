package com.github.zimmerlab.gtfcompare.utils;

import com.github.kleinsamuel.gtfutils.GtfFile;
import com.github.kleinsamuel.gtfutils.feature.GeneFeature;
import com.github.zimmerlab.gtfcompare.model.GenePair;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import java.util.*;
import java.util.function.Function;

public class OverlappingGenes {

    public static Map<GeneFeature, GeneFeature> map(GtfFile targetGtfFile, GtfFile queryGtfFile) {
        var overlaps = findOverlaps(targetGtfFile, queryGtfFile);

        var locusPairs = findLoci(overlaps);

        var finalMapping = new LinkedHashMap<GeneFeature, GeneFeature>();
        var unresolvedPairs = new ArrayList<GenePair>();
        for (var locus : locusPairs) {
            Map<GeneFeature, GeneFeature> exact = findExactMatches(locus);
            finalMapping.putAll(exact);

            var unresolved = locus.stream().filter(p -> !exact.containsKey(p.getTargetGene())).toList();
            unresolvedPairs.addAll(unresolved);
        }

        return finalMapping;
    }


    public static List<GenePair> findOverlaps(GtfFile targetGtfFile, GtfFile queryGtfFile) {
        var targetTrees = new HashMap<String, IntervalTreeMap<GeneFeature>>();
        var queryTrees = new HashMap<String, IntervalTreeMap<GeneFeature>>();

        Function<GeneFeature, String> keyOf = g -> g.getBaseData().getContig();

        extractGenes(targetGtfFile, targetTrees, keyOf);
        extractGenes(queryGtfFile, queryTrees, keyOf);

        var pairs = new ArrayList<GenePair>();
        for (var e : targetTrees.entrySet()) {
            String key = e.getKey();
            var tTree = e.getValue();
            var qTree = queryTrees.get(key);

            if (qTree == null) {
                continue;
            }

            for (GeneFeature t : tTree.values()) {
                var td = t.getBaseData();
                var tiv = new Interval(td.getContig(), td.getStart(), td.getEnd());

                for (GeneFeature q : qTree.getOverlapping(tiv)) {
                    pairs.add(new GenePair(t, q));
                }
            }
        }

        return pairs;
    }

    private static void extractGenes(GtfFile gtfFile, HashMap<String, IntervalTreeMap<GeneFeature>> trees, Function<GeneFeature, String> keyOf) {
        for (var id : gtfFile.getAllGeneFeatureIds()) {
            var g = gtfFile.getGeneFeature(id);
            var d = g.getBaseData();
            String key = keyOf.apply(g);
            trees.computeIfAbsent(key, k -> new IntervalTreeMap<>()).put(new Interval(d.getContig(), d.getStart(), d.getEnd()), g);
        }
    }

    private static class DisjointSet<T> {
        private final Map<T, T> parent = new HashMap<>();

        T find(T x) {
            parent.putIfAbsent(x, x);
            if (!parent.get(x).equals(x)) {
                parent.put(x, find(parent.get(x)));
            }
            return parent.get(x);
        }

        void union(T a, T b) {
            T ra = find(a), rb = find(b);
            if (!ra.equals(rb)) parent.put(ra, rb);
        }
    }

    public static List<List<GenePair>> findLoci(List<GenePair> pairs) {
        var ds = new DisjointSet<GeneFeature>();
        for (var p : pairs) {
            var t = p.getTargetGene();
            var q = p.getQueryGene();
            ds.union(t, q);
        }

        var rootToPairs = new HashMap<GeneFeature, List<GenePair>>();
        for (var p : pairs) {
            GeneFeature root = ds.find(p.getTargetGene());
            rootToPairs.computeIfAbsent(root, k -> new ArrayList<>()).add(p);
        }

        return new ArrayList<>(rootToPairs.values());
    }

    private static Map<GeneFeature, GeneFeature> findExactMatches(List<GenePair> locus) {
        var exact = new HashMap<GeneFeature, GeneFeature>();
        for (var p : locus) {
            var t = p.getTargetGene().getBaseData();
            var q = p.getQueryGene().getBaseData();
            if (t.getStart() == q.getStart() && t.getEnd() == q.getEnd()) {
                exact.put(p.getTargetGene(), p.getQueryGene());
            }
        }
        return exact;
    }
}
