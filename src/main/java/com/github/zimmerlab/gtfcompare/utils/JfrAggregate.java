package com.github.zimmerlab.gtfcompare.utils;

import jdk.jfr.consumer.RecordedEvent;
import jdk.jfr.consumer.RecordedFrame;
import jdk.jfr.consumer.RecordedMethod;
import jdk.jfr.consumer.RecordingFile;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.time.Duration;
import java.time.Instant;
import java.util.*;

public class JfrAggregate {
    private static final Logger logger = LoggerFactory.getLogger(JfrAggregate.class);

    private static final String EV_EXEC = "jdk.ExecutionSample";
    private static final String EV_ALLOC_SAMPLE = "jdk.ObjectAllocationSample";

    private static final int TOP = 100;
    private static long totalCpuSamples = 0L;
    private static long totalAllocEvents = 0L;

    private static class Counters {
        long cpuSamples = 0L;
        long allocBytes = 0L;

        void addCpu(long v) {
            cpuSamples += v;
        }

        void addAlloc(long b) {
            allocBytes += b;
        }
    }

    public static void benchmark(String jfrFile, String outDir, String prefix) throws Exception {
        var byMethod = new HashMap<String, Counters>();
        long totalAllocatedBytes = 0L;
        var totalDuration = Duration.ZERO;

        boolean sawAllocSample = false;
        boolean sawAllocTlabs = false;

        var p = Paths.get(jfrFile);
        if (!Files.isRegularFile(p)) {
            System.err.println("Not a file: " + jfrFile);
            return;
        }

        try (var rf = new RecordingFile(p)) {
            Instant minStart = null;
            Instant maxEnd = null;

            while (rf.hasMoreEvents()) {
                var e = rf.readEvent();

                Instant t = e.getStartTime();
                if (t != null) {
                    if (minStart == null || t.isBefore(minStart)) {
                        minStart = t;
                    }
                    if (maxEnd == null || t.isAfter(maxEnd)) {
                        maxEnd = t;
                    }
                }

                String name = e.getEventType().getName();

                switch (name) {
                    case EV_EXEC -> {
                        totalCpuSamples++;
                        handleExecSample(e, byMethod);
                    }
                    case EV_ALLOC_SAMPLE-> {
                        totalAllocEvents++;
                        sawAllocSample = true;
                        long size = getAllocationBytes(e);
                        totalAllocatedBytes += size;
                        attributeAllocToMethod(e, size, byMethod);
                    }
                    default -> {}
                }
            }

            if (minStart != null && maxEnd != null && maxEnd.isAfter(minStart)) {
                totalDuration = totalDuration.plus(Duration.between(minStart, maxEnd));
            }

        } catch (IOException ioe) {
            System.err.println("Failed to read " + jfrFile + ": " + ioe.getMessage());
        }



        var topCpu = byMethod.entrySet().stream().sorted(Comparator.comparingLong((Map.Entry<String, Counters> e) -> e.getValue().cpuSamples).reversed()).limit(TOP).toList();
        var topAlloc = byMethod.entrySet().stream().sorted(Comparator.comparingLong((Map.Entry<String, Counters> e) -> e.getValue().allocBytes).reversed()).limit(TOP).toList();

        try {
            writeTsvReports(
                    Paths.get(outDir),
                    totalDuration, totalAllocatedBytes,
                    sawAllocSample, sawAllocTlabs,
                    totalCpuSamples, totalAllocEvents,
                    topCpu, topAlloc, byMethod, prefix
            );
        } catch (IOException ex) {
            logger.error("Failed to write TSV reports: ", ex);
        }
    }

    private static void handleExecSample(RecordedEvent e, Map<String, Counters> byMethod) {
        String key = topUserOrJavaFrameKey(e);
        if (key == null) return;
        byMethod.computeIfAbsent(key, k -> new Counters()).addCpu(1L);
    }

    private static void attributeAllocToMethod(RecordedEvent e, long size, Map<String, Counters> byMethod) {
        String key = topUserOrJavaFrameKey(e);
        if (key == null) return;
        byMethod.computeIfAbsent(key, k -> new Counters()).addAlloc(size);
    }

    private static String methodKey(RecordedMethod m, RecordedFrame f) {
        var owner = (m.getType() != null) ? m.getType().getName() : "<unknown>";
        var name = m.getName();
        var sig = toJavaSignature(m.getDescriptor());
        //Integer line = f.getLineNumber();
        return owner + "." + name + sig; // + (line != null ? (":" + line) : "");
    }

    // Choose the first frame in the stack that belongs to "user code" (by prefix).
// If none is found, fall back to the first Java frame.
    private static String topUserOrJavaFrameKey(RecordedEvent e) {
        var st = e.getStackTrace();
        if (st == null || st.getFrames() == null) return null;

        String fallback = null;

        for (var f : st.getFrames()) {
            var m = f.getMethod();
            if (m == null || m.getType() == null) continue;
            String typeName = m.getType().getName();

            // remember first Java frame as fallback
            if (fallback == null) {
                fallback = methodKey(m, f);
            }

            if (startsWithAny(typeName, EXCLUDE_PACKAGES)) continue;

            // prefer user packages
            if (startsWithAny(typeName, INCLUDE_PACKAGES)) {
                return methodKey(m, f);
            }
        }
        return fallback;
    }


    private static void printTable(List<Map.Entry<String, Counters>> rows, boolean showCpu, boolean showAlloc) {
        System.out.printf("%-8s  %-12s  %s%n", showCpu ? "CPU" : "", showAlloc ? "Alloc(Bytes)" : "", "Method");
        for (var r : rows) {
            String left = showCpu ? String.format("%-8d", r.getValue().cpuSamples) : "";
            String mid = showAlloc ? String.format("%-12d", r.getValue().allocBytes) : "";
            System.out.printf("%s  %s  %s%n", left, mid, r.getKey());
        }
    }

    private static String human(Duration d) {
        long s = d.getSeconds();
        long h = s / 3600;
        s %= 3600;
        long m = s / 60;
        s %= 60;
        long ms = d.minusSeconds(d.getSeconds()).toMillis();
        return String.format("%02dh:%02dm:%02ds.%03d", h, m, s, ms);
    }

    private static String csvEscape(String s) {
        if (s.contains(",") || s.contains("\"")) {
            return "\"" + s.replace("\"", "\"\"") + "\"";
        }
        return s;
    }

    private static long getAllocationBytes(RecordedEvent e) {
        String[] candidates = {"allocationSize", "size", "weight"};
        for (String f : candidates) {
            if (e.getEventType().getField(f) != null) {
                Object v = e.getValue(f);
                if (v instanceof Number num) {
                    long bytes = (long) Math.max(0, num.doubleValue());
                    return bytes;
                }
            }
        }
        return 0L;
    }

    private static final List<String> INCLUDE_PACKAGES = List.of("com.github.zimmerlab.", "com.github.kleinsamuel.");

    private static final List<String> EXCLUDE_PACKAGES = List.of(
            "java.", "jdk.", "sun.", "javax.",
            "org.springframework.", "org.apache.", "com.fasterxml.",
            "kotlin.", "scala."
    );

    private static boolean startsWithAny(String s, List<String> prefixes) {
        for (String p : prefixes) {
            if (s.startsWith(p)) return true;
        }
        return false;
    }

    private static String toJavaSignature(String desc) {
        if (desc == null || desc.isEmpty()) return "()";
        int endParams = desc.indexOf(')');
        if (desc.charAt(0) != '(' || endParams < 0) return desc;

        var params = desc.substring(1, endParams);
        var ret = desc.substring(endParams + 1);

        var ptypes = new ArrayList<String>();
        for (int i = 0; i < params.length(); ) {
            int[] r = parseType(params, i);
            ptypes.add(jvmToJavaType(params.substring(i, r[0])));
            i = r[0];
        }
        return "(" + String.join(", ", ptypes) + ")" + jvmToJavaType(ret);
    }

    private static int[] parseType(String s, int i) {
        int dims = 0;
        while (i < s.length() && s.charAt(i) == '[') {
            dims++;
            i++;
        }
        char c = s.charAt(i);
        if (c == 'L') {
            int semi = s.indexOf(';', i);
            i = semi + 1;
        } else {
            i++;
        }
        return new int[]{i, dims};
    }

    private static String jvmToJavaType(String t) {
        int dims = 0;
        while (dims < t.length() && t.charAt(dims) == '[') dims++;
        String core = t.substring(dims);
        String base = switch (core.charAt(0)) {
            case 'V' -> "void";
            case 'Z' -> "boolean";
            case 'B' -> "byte";
            case 'C' -> "char";
            case 'S' -> "short";
            case 'I' -> "int";
            case 'J' -> "long";
            case 'F' -> "float";
            case 'D' -> "double";
            case 'L' -> core.substring(1, core.length() - 1).replace('/', '.');
            default -> core;
        };
        return base + "[]".repeat(dims);
    }

    private static void writeTsvReports(
            java.nio.file.Path outDir,
            java.time.Duration totalDuration,
            long totalAllocatedBytes,
            boolean sawAllocSample,
            boolean sawAllocTlabs,
            long totalCpuSamples,
            long totalAllocEvents,
            java.util.List<java.util.Map.Entry<String, Counters>> topCpu,
            java.util.List<java.util.Map.Entry<String, Counters>> topAlloc,
            java.util.Map<String, Counters> byMethod,
            String prefix
    ) throws java.io.IOException {
        java.nio.file.Files.createDirectories(outDir);

        var summaryName = "summary.tsv";
        var top_cpu_name = "top_cpu.tsv";
        var top_alloc_name = "top_alloc.tsv";
        var methods_name = "methods.tsv";

        if(prefix != null){
            summaryName = prefix + "_" + summaryName;
            top_cpu_name = prefix + "_" + top_cpu_name;
            top_alloc_name = prefix + "_" + top_alloc_name;
            methods_name = prefix + "_" + methods_name;
        }
        try (var w = java.nio.file.Files.newBufferedWriter(outDir.resolve(summaryName))) {
            // header
            w.write("totalRecordingsDurationMillis\ttotalAllocatedBytes\tallocationSources\ttotalCpuSamples\ttotalAllocationEvents\n");
            // single row
            String allocSources = (sawAllocSample ? "ObjectAllocationSample " : "")
                    + (sawAllocTlabs ? "ObjectAllocationIn/OutsideTLAB"
                    : (sawAllocSample ? "" : "(none)"));
            w.write(totalDuration.toMillis() + "\t"
                    + totalAllocatedBytes + "\t"
                    + tsvEscape(allocSources) + "\t"
                    + totalCpuSamples + "\t"
                    + totalAllocEvents + "\n");
        }

        try (var w = java.nio.file.Files.newBufferedWriter(outDir.resolve(top_cpu_name))) {
            w.write("rank\tmethod\tcpuSamples\n");
            for (int i = 0; i < topCpu.size(); i++) {
                var e = topCpu.get(i);
                w.write((i + 1) + "\t" + tsvEscape(e.getKey()) + "\t" + e.getValue().cpuSamples + "\n");
            }
        }

        try (var w = java.nio.file.Files.newBufferedWriter(outDir.resolve(top_alloc_name))) {
            w.write("rank\tmethod\tallocBytes\n");
            for (int i = 0; i < topAlloc.size(); i++) {
                var e = topAlloc.get(i);
                w.write((i + 1) + "\t" + tsvEscape(e.getKey()) + "\t" + e.getValue().allocBytes + "\n");
            }
        }

        try (var w = java.nio.file.Files.newBufferedWriter(outDir.resolve(methods_name))) {
            w.write("method\tcpuSamples\tallocBytes\n");
            byMethod.entrySet().stream()
                    .sorted(java.util.Comparator.comparingLong(
                            (java.util.Map.Entry<String, Counters> e) -> e.getValue().cpuSamples + e.getValue().allocBytes
                    ).reversed())
                    .forEach(e -> {
                        try {
                            w.write(tsvEscape(e.getKey()) + "\t" + e.getValue().cpuSamples + "\t" + e.getValue().allocBytes + "\n");
                        } catch (java.io.IOException ex) {
                            throw new RuntimeException(ex);
                        }
                    });
        }
    }

    private static String tsvEscape(String s) {
        if (s == null) return "";
        String r = s.replace('\t', ' ');
        r = r.replace('\r', ' ').replace('\n', ' ');
        return r;
    }
}