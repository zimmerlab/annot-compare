package com.github.zimmerlab.gtfcompare.analysis.clique;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.stream.Collectors;

public class Reporter {
    private final List<AnalysisResult> results = new ArrayList<>();

    public void add(AnalysisResult r) { results.add(r); }

    public void printTable(PrintStream out) {
        // Fixed-width layout; adjust widths to your data
        String headerFmt = "%-15s %8s %8s %9s %12s %12s %12s %18s %18s%n";
        String rowFmt    = "%-15s %8d %8d %9d %12d %12d %12d %18d %18d%n";
        out.printf(headerFmt, "Contig", "Target", "Query", "Clusters", "ExactMatch",
                "Unm.Tgt", "Unm.Qry", "TgtNoExact(inC)", "QryNoExact(inC)");
        out.println("-----------------------------------------------------------------------------------------------"
                + "-------------------------------------------------------");
        for (var r : results) {
            out.printf(rowFmt, r.contig, r.targetGenes, r.queryGenes, r.clusters, r.exactMatches,
                    r.unmatchedTargetGenes, r.unmatchedQueryGenes,
                    r.inClustersNoExactTarget, r.inClustersNoExactQuery);
        }
        var tot = summarize();
        out.println();
        out.printf(rowFmt, "TOTAL",
                tot.targetGenes, tot.queryGenes, tot.clusters, tot.exactMatches,
                tot.unmatchedTargetGenes, tot.unmatchedQueryGenes,
                tot.inClustersNoExactTarget, tot.inClustersNoExactQuery);
    }

    private AnalysisResult summarize() {
        return new AnalysisResult(
                "TOTAL",
                results.stream().mapToInt(r -> r.targetGenes).sum(),
                results.stream().mapToInt(r -> r.queryGenes).sum(),
                results.stream().mapToInt(r -> r.clusters).sum(),
                results.stream().mapToInt(r -> r.exactMatches).sum(),
                results.stream().mapToInt(r -> r.unmatchedTargetGenes).sum(),
                results.stream().mapToInt(r -> r.unmatchedQueryGenes).sum(),
                results.stream().mapToInt(r -> r.inClustersNoExactTarget).sum(),
                results.stream().mapToInt(r -> r.inClustersNoExactQuery).sum()
        );
    }

    public void writeCsv(Path path) throws IOException {
        try (var w = Files.newBufferedWriter(path, StandardCharsets.UTF_8)) {
            w.write("contig,targetGenes,queryGenes,clusters,exactMatches,unmatchedTarget,unmatchedQuery,"
                    + "inClustersNoExactTarget,inClustersNoExactQuery\n");
            for (var r : results) {
                w.write(String.format(Locale.ROOT,
                        "%s,%d,%d,%d,%d,%d,%d,%d,%d%n",
                        r.contig, r.targetGenes, r.queryGenes, r.clusters, r.exactMatches,
                        r.unmatchedTargetGenes, r.unmatchedQueryGenes,
                        r.inClustersNoExactTarget, r.inClustersNoExactQuery));
            }
        }
    }
}
