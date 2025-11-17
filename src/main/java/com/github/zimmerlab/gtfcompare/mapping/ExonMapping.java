package com.github.zimmerlab.gtfcompare.mapping;

import com.github.kleinsamuel.gtfutils.GtfConfig;
import com.github.kleinsamuel.gtfutils.feature.GtfFeature;
import com.github.kleinsamuel.gtfutils.feature.TranscriptFeature;
import com.github.zimmerlab.gtfcompare.model.FeaturePair;
import com.github.zimmerlab.gtfcompare.utils.Constants;

import java.util.*;
import java.util.function.BiPredicate;
import java.util.function.IntSupplier;
import java.util.function.ToIntFunction;

public class ExonMapping {
    public static List<GtfFeature> sortedExons(TranscriptFeature transcriptFeature){
        return sortedExons(transcriptFeature, true);
    }

    public static List<GtfFeature> sortedExons(TranscriptFeature tf, boolean sortByStrand) {
        boolean fwd = tf.getBaseData().isForwardStrand();
        Comparator<GtfFeature> byStart = Comparator.comparingInt(f -> f.getBaseData().getStart());
        var list = tf.getFeatures().stream().filter(f -> "exon".equals(f.getBaseData().getType())).sorted(byStart).toList();

        if (sortByStrand) {
            return fwd ? list : new ArrayList<>(list) {{
                Collections.reverse(this);
            }};
        }

        return list;
    }

    private static List<Integer> gapProfile(List<GtfFeature> exons) {
        var gaps = new ArrayList<Integer>(Math.max(0, exons.size() - 1));
        for (int i = 0; i < exons.size() - 1; i++) {
            int g = exons.get(i + 1).getBaseData().getStart() - exons.get(i).getBaseData().getEnd() - 1;
            gaps.add(g);
        }
        return gaps;
    }

    private static int[][] nwBacktrace(int[] A, int[] B, int gapPenalty, int cap) {
        int n = A.length, m = B.length;
        var dp = new int[n + 1][m + 1];
        var bt = new int[n + 1][m + 1]; // 0=diag,1=up,2=left

        for (int i = 1; i <= n; i++) {
            dp[i][0] = dp[i - 1][0] - gapPenalty;
            bt[i][0] = 1;
        }
        for (int j = 1; j <= m; j++) {
            dp[0][j] = dp[0][j - 1] - gapPenalty;
            bt[0][j] = 2;
        }

        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                int diff = Math.abs(A[i - 1] - B[j - 1]);
                int match = dp[i - 1][j - 1] - Math.min(diff, cap); // diag
                int delA = dp[i - 1][j] - gapPenalty; // up
                int delB = dp[i][j - 1] - gapPenalty; // left

                if (match >= delA && match >= delB) {
                    bt[i][j] = 0;
                    dp[i][j] = match;
                } else if (delA >= delB) {
                    bt[i][j] = 1;
                    dp[i][j] = delA;
                } else {
                    bt[i][j] = 2;
                    dp[i][j] = delB;
                }
            }
        }
        return bt;
    }

    private static List<int[]> recoverGapMatches(int[][] bt, int n, int m) {
        var pairs = new ArrayList<int[]>();
        int i = n, j = m;
        while (i > 0 || j > 0) {
            if (i > 0 && j > 0 && bt[i][j] == 0) {
                pairs.add(new int[]{i - 1, j - 1});
                i--;
                j--;
            } else if (i > 0 && (j == 0 || bt[i][j] == 1)) {
                i--;
            } else {
                j--;
            }
        }
        Collections.reverse(pairs);
        return pairs;
    }

    private static int endpointManhattan(GtfFeature a, GtfFeature b) {
        return Math.abs(a.getBaseData().getStart() - b.getBaseData().getStart()) + Math.abs(a.getBaseData().getEnd() - b.getBaseData().getEnd());
    }

    public static List<FeaturePair<GtfFeature>> pairExonsByGapAlignment(TranscriptFeature t, TranscriptFeature q, int gapPenalty, int capDelta, int lenCapBp, double lenCapFrac) {
        var pairs = new ArrayList<FeaturePair<GtfFeature>>();
        var targetExons = sortedExons(t);
        var queryExons = sortedExons(q);

        if (targetExons.isEmpty() && queryExons.isEmpty()) return pairs;
        if (targetExons.size() == 1 && queryExons.size() == 1) {
            pairs.add(new FeaturePair<GtfFeature>(targetExons.getFirst(), queryExons.getFirst()));
            return pairs;
        }

        var A = gapProfile(targetExons).stream().mapToInt(x -> x).toArray();
        var B = gapProfile(queryExons).stream().mapToInt(x -> x).toArray();

        ToIntFunction<GtfFeature> len = f -> f.getBaseData().getEnd() - f.getBaseData().getStart() + 1;

        BiPredicate<GtfFeature, GtfFeature> passLen = (te, qe) -> {
            int lt = len.applyAsInt(te), lq = len.applyAsInt(qe);
            int abs = Math.abs(lt - lq);
            double rel = (double) abs / Math.max(lt, lq);
            return abs <= lenCapBp || rel <= lenCapFrac;
        };

        if (A.length == 0 || B.length == 0) {
            var usedQ = new HashSet<GtfFeature>();
            for (var te : targetExons) {
                GtfFeature best = null;
                int bestLenDiff = Integer.MAX_VALUE;
                int bestD = Integer.MAX_VALUE;
                for (var qe : queryExons) {
                    if (usedQ.contains(qe)) continue;
                    if (!passLen.test(te, qe)) continue;
                    int lenDiff = Math.abs(len.applyAsInt(te) - len.applyAsInt(qe));
                    int d = endpointManhattan(te, qe);
                    if (lenDiff < bestLenDiff || (lenDiff == bestLenDiff && d < bestD)) {
                        bestLenDiff = lenDiff;
                        bestD = d;
                        best = qe;
                    }
                }
                pairs.add(new FeaturePair<>(te, best));
                if (best != null) usedQ.add(best);
            }
            for (var qe : queryExons)
                if (pairs.stream().noneMatch(p -> qe.equals(p.getQuery()))) pairs.add(new FeaturePair<>(null, qe));
            return pairs;
        }

        var bt = nwBacktrace(A, B, gapPenalty, capDelta);
        var gapMatches = recoverGapMatches(bt, A.length, B.length);

        var usedT = new boolean[targetExons.size()];
        var usedQ = new boolean[queryExons.size()];

        var deltas = new ArrayList<Integer>();
        IntSupplier deltaHat = () -> {
            if (deltas.size() < 2) return 0;
            var copy = new ArrayList<>(deltas);
            copy.sort(Integer::compare);
            int mid = copy.size() / 2;
            return (copy.size() % 2 == 1) ? copy.get(mid) : (copy.get(mid - 1) + copy.get(mid)) / 2;
        };

        for (var g : gapMatches) {
            var k = g[0];
            var l = g[1];
            int[] candT = {k, k + 1};
            int[] candQ = {l, l + 1};

            record Cand(int ct, int cq, int lenDiff, int rank) {
            }
            var cands = new ArrayList<Cand>();

            int delta = deltaHat.getAsInt();

            for (int ct : candT) {
                for (int cq : candQ) {
                    if (ct < 0 || cq < 0 || ct >= targetExons.size() || cq >= queryExons.size()) continue;
                    if (usedT[ct] || usedQ[cq]) continue;
                    var te = targetExons.get(ct);
                    var qe = queryExons.get(cq);
                    if (!passLen.test(te, qe)) continue;

                    int lenDiff = Math.abs(len.applyAsInt(te) - len.applyAsInt(qe));

                    int rank = Math.abs((te.getBaseData().getStart() - qe.getBaseData().getStart()) - delta) + Math.abs((te.getBaseData().getEnd() - qe.getBaseData().getEnd()) - delta);

                    cands.add(new Cand(ct, cq, lenDiff, rank));
                }
            }

            cands.sort((x, y) -> {
                if (x.lenDiff != y.lenDiff) return Integer.compare(x.lenDiff, y.lenDiff); // 1) min lenDiff
                return Integer.compare(x.rank, y.rank);                                    // 2) min rank
            });

            if (!cands.isEmpty()) {
                var best = cands.get(0);
                var te = targetExons.get(best.ct);
                var qe = queryExons.get(best.cq);
                pairs.add(new FeaturePair<>(te, qe));
                usedT[best.ct] = usedQ[best.cq] = true;

                deltas.add(te.getBaseData().getStart() - qe.getBaseData().getStart());
            }
        }

// right anchor
        int iR = targetExons.size() - 1, jR = queryExons.size() - 1;
        while (iR >= 0 && usedT[iR]) iR--;
        while (jR >= 0 && usedQ[jR]) jR--;
        if (iR >= 0 && jR >= 0) {
            var te = targetExons.get(iR);
            var qe = queryExons.get(jR);
            int lt = len.applyAsInt(te);
            int lq = len.applyAsInt(qe);
            int lenDiff = Math.abs(lt - lq);
            double rel = (double) lenDiff / Math.max(lt, lq);
            if (lenDiff <= lenCapBp || rel <= lenCapFrac) {
                pairs.add(new FeaturePair<>(te, qe));
                usedT[iR] = usedQ[jR] = true;
            }
        }

// greedy for the rest
        int delta = 0;
        if (!pairs.isEmpty()) {
            var ds = new ArrayList<Integer>();
            for (var p : pairs) {
                if (p.getTarget() != null && p.getQuery() != null) {
                    ds.add(p.getTarget().getBaseData().getStart() - p.getQuery().getBaseData().getStart());
                }
            }
            if (ds.size() >= 2) {
                ds.sort(Integer::compare);
                int mid = ds.size() / 2;
                delta = (ds.size() % 2 == 1) ? ds.get(mid) : (ds.get(mid - 1) + ds.get(mid)) / 2;
            }
        }
        for (int it = 0; it < targetExons.size(); it++) {
            if (usedT[it]) continue;
            var te = targetExons.get(it);
            GtfFeature best = null;
            int bestLenDiff = Integer.MAX_VALUE;
            int bestRank = Integer.MAX_VALUE;
            for (int iq = 0; iq < queryExons.size(); iq++) {
                if (usedQ[iq]) continue;
                var qe = queryExons.get(iq);
                int lt = len.applyAsInt(te);
                int lq = len.applyAsInt(qe);

                int lenDiff = Math.abs(lt - lq);
                double rel = (double) lenDiff / Math.max(lt, lq);
                if (!(lenDiff <= lenCapBp || rel <= lenCapFrac)) continue;

                int rank = Math.abs((te.getBaseData().getStart() - qe.getBaseData().getStart()) - delta) + Math.abs((te.getBaseData().getEnd() - qe.getBaseData().getEnd()) - delta);
                if (lenDiff < bestLenDiff || (lenDiff == bestLenDiff && rank < bestRank)) {
                    bestLenDiff = lenDiff;
                    bestRank = rank;
                    best = qe;
                }
            }
            if (best != null) {
                pairs.add(new FeaturePair<>(te, best));
                int idx = queryExons.indexOf(best);
                usedT[it] = true;
                usedQ[idx] = true;
            }
        }

        for (int i = 0; i < targetExons.size(); i++)
            if (!usedT[i]) pairs.add(new FeaturePair<>(targetExons.get(i), null));
        for (int j = 0; j < queryExons.size(); j++)
            if (!usedQ[j]) pairs.add(new FeaturePair<>(null, queryExons.get(j)));

        return pairs;
    }


    public static boolean within(int s, int e, int xs, int xe, int pad) {
        return s >= xs - pad && e <= xe + pad;
    }

    public static List<GtfFeature> featuresOfType(TranscriptFeature tf, String type) {
        return tf.getFeatures().stream().filter(f -> type.equals(GtfConfig.getDefault(f.getBaseData().getType()))).sorted(Comparator.comparingInt(f -> f.getBaseData().getStart())).toList();
    }

    public static List<FeaturePair<GtfFeature>> mapFeaturesWithinExonPairs(TranscriptFeature t, TranscriptFeature q, List<FeaturePair<GtfFeature>> exonPairs, String featureType, int padBp) {
        var res = new ArrayList<FeaturePair<GtfFeature>>();
        var usedQ = new HashSet<GtfFeature>();

        var tFeat = featuresOfType(t, featureType);
        var qFeat = featuresOfType(q, featureType);

        var qByExon = new HashMap<GtfFeature, List<GtfFeature>>();
        for (var fp : exonPairs) {
            var qe = fp.getQuery();
            if (qe != null) qByExon.put(qe, new ArrayList<>());
        }
        for (var qf : qFeat) {
            for (var fp : exonPairs) {
                var qe = fp.getQuery();
                if (qe == null) continue;
                var qS = qf.getBaseData().getStart();
                var qE = qf.getBaseData().getEnd();
                var qeS = fp.getQuery().getBaseData().getStart();
                var qeE = fp.getQuery().getBaseData().getEnd();
                if (within(qS, qE, qeS, qeE, padBp)) {
                    qByExon.computeIfAbsent(qe, k -> new ArrayList<>()).add(qf);
                    break;
                }
            }
        }

        for (var fpExon : exonPairs) {
            var targetExon = fpExon.getTarget();
            var queryExon = fpExon.getQuery();
            var tSub = new ArrayList<GtfFeature>();
            if (targetExon != null) {
                int teS = targetExon.getBaseData().getStart(), teE = targetExon.getBaseData().getEnd();
                for (var tf : tFeat) {
                    int s = tf.getBaseData().getStart(), e = tf.getBaseData().getEnd();
                    if (within(s, e, teS, teE, padBp)) tSub.add(tf);
                }
            }

            var candidates = (queryExon != null) ? qByExon.getOrDefault(queryExon, List.of()) : List.<GtfFeature>of();

            for (var tf : tSub) {
                GtfFeature best = null;
                int bestDist = Integer.MAX_VALUE;
                for (var qf : candidates) {
                    if (usedQ.contains(qf)) continue;
                    int d = Math.abs(tf.getBaseData().getStart() - qf.getBaseData().getStart()) + Math.abs(tf.getBaseData().getEnd() - qf.getBaseData().getEnd());
                    if (d < bestDist) {
                        bestDist = d;
                        best = qf;
                    }
                }
                res.add(new FeaturePair<>(tf, best));
                if (best != null) usedQ.add(best);
            }
        }

        for (var qf : qFeat) {
            if (!usedQ.contains(qf) && res.stream().noneMatch(p -> qf.equals(p.getQuery()))) {
                res.add(new FeaturePair<>(null, qf));
            }
        }
        return res;
    }

    public static List<FeaturePair<GtfFeature>> mapIntronsByExonPairs(
            TranscriptFeature targetTranscript,
            TranscriptFeature queryTranscript,
            List<FeaturePair<GtfFeature>> exonPairs
    ) {
        var result = new ArrayList<FeaturePair<GtfFeature>>();
        var usedT = new HashSet<GtfFeature>();
        var usedQ = new HashSet<GtfFeature>();

        var tIntrons = featuresOfType(targetTranscript, Constants.INTRON);
        var qIntrons = featuresOfType(queryTranscript, Constants.INTRON);

        var tExonsOrdered = sortedExons(targetTranscript, true);
        var qExonsOrdered = sortedExons(queryTranscript, true);

        var tIndex = new HashMap<GtfFeature, Integer>();
        for (int i = 0; i < tExonsOrdered.size(); i++) tIndex.put(tExonsOrdered.get(i), i);

        var qIndex = new HashMap<GtfFeature, Integer>();
        for (int i = 0; i < qExonsOrdered.size(); i++) qIndex.put(qExonsOrdered.get(i), i);

        var t2q = new HashMap<GtfFeature, GtfFeature>();
        for (var fp : exonPairs) {
            var te = fp.getTarget();
            var qe = fp.getQuery();
            if (te != null && qe != null) {
                t2q.put(te, qe);
            }
        }

        for (int i = 0; i + 1 < tExonsOrdered.size(); i++) {
            var targetExonLeft = tExonsOrdered.get(i);
            var targetExonRight = tExonsOrdered.get(i + 1);

            var queryExonLeft = t2q.get(targetExonLeft);
            var queryExonRight = t2q.get(targetExonRight);
            if (queryExonLeft == null || queryExonRight == null) {
                continue;
            }

            var queryIndexLeft = qIndex.get(queryExonLeft);
            var queryIndexRight = qIndex.get(queryExonRight);
            if (queryIndexLeft == null || queryIndexRight == null) continue;

            if (queryIndexRight != queryIndexLeft + 1) {
                continue;
            }

            var tIntron = findIntronBetween(targetExonLeft, targetExonRight, tIntrons);
            var qIntron = findIntronBetween(queryExonLeft, queryExonRight, qIntrons);

            if (tIntron != null || qIntron != null) {
                result.add(new FeaturePair<>(tIntron, qIntron));
                if (tIntron != null) usedT.add(tIntron);
                if (qIntron != null) usedQ.add(qIntron);
            }
        }

        for (var ti : tIntrons) {
            if (!usedT.contains(ti)) result.add(new FeaturePair<>(ti, null));
        }
        for (var qi : qIntrons) {
            if (!usedQ.contains(qi)) result.add(new FeaturePair<>(null, qi));
        }

        return result;
    }

    private static GtfFeature findIntronBetween(
            GtfFeature eLeft,
            GtfFeature eRight,
            List<GtfFeature> introns
    ) {
        int exonLeftEnd = eLeft.getBaseData().getEnd() + 1;
        int exonRightStart = eRight.getBaseData().getStart() - 1;
        if (exonRightStart < exonLeftEnd) return null;

        GtfFeature best = null;
        int bestOverlap = -1;

        for (var intron : introns) {
            int start = intron.getBaseData().getStart();
            int end = intron.getBaseData().getEnd();

            if (!within(start, end, exonLeftEnd, exonRightStart, 0)) continue;

            int overlap = Math.min(end, exonRightStart) - Math.max(start, exonLeftEnd) + 1;
            if (overlap > bestOverlap) {
                bestOverlap = overlap;
                best = intron;
            }
        }
        return best;
    }

    // DEPRECATED

    private List<FeaturePair> pairByPosition(List<GtfFeature> targets, List<GtfFeature> queries) {

        Comparator<GtfFeature> byStartThenEnd = Comparator.comparingInt((GtfFeature f) -> f.getBaseData().getStart()).thenComparingInt(f -> f.getBaseData().getEnd());

        var a = new ArrayList<>(targets);
        var b = new ArrayList<>(queries);
        a.sort(byStartThenEnd);
        b.sort(byStartThenEnd);

        var aSize = a.size();
        var bSize = b.size();
        var n = Math.max(aSize, bSize);
        var pairs = new ArrayList<FeaturePair>(n);
        for (int i = 0; i < n; i++) {
            GtfFeature ta = i < aSize ? a.get(i) : null;
            GtfFeature qb = i < bSize ? b.get(i) : null;
            pairs.add(new FeaturePair(ta, qb));
        }
        return pairs;
    }

    private List<FeaturePair> pairByOverlap(List<GtfFeature> targets, List<GtfFeature> queries) {
        var pairs = new ArrayList<FeaturePair>();
        for (GtfFeature ta : targets) {
            if (ta.getBaseData().getAttributes("transcript_id") != null && ta.getBaseData().getAttributes("transcript_id").getFirst().equals("ENST00000412513")) {
                var a = 2;
            }
            var targetStart = ta.getBaseData().getStart();
            var targetEnd = ta.getBaseData().getEnd();

            var found = false;
            for (GtfFeature qb : queries) {
                var queryStart = qb.getBaseData().getStart();
                int queryEnd = qb.getBaseData().getEnd();

                if (queryStart <= targetEnd && queryEnd >= targetStart) {
                    pairs.add(new FeaturePair(ta, qb));
                    found = true;
                }
            }

            if (!found) {
                pairs.add(new FeaturePair(ta, null));
            }
        }

        for (GtfFeature qb : queries) {
            var matched = pairs.stream().anyMatch(p -> p.getQuery() != null && p.getQuery().equals(qb));
            if (!matched) {
                pairs.add(new FeaturePair(null, qb));
            }
        }

        return pairs;
    }
}
