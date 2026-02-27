// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
// Location: MassSpectrometry/Dia/Scoring/DiaFdrEngine.cs

using System;
using System.Buffers;
using System.Collections.Generic;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Target-decoy FDR estimation engine for DIA precursor identifications.
    /// 
    /// Implements the iterative semi-supervised learning approach used by DIA-NN:
    ///   Iteration 1: Train on pseudo-labels — top 40% targets (positive) vs ALL decoys (negative).
    ///   Iteration 2+: Retrain using q &lt; 0.01 targets (positive) vs subsampled decoys (negative).
    ///   Repeat until convergence (weight change &lt; 1%, or ID count stable, or max 5 iterations).
    /// 
    /// Q-value calculation:
    ///   - Sort by classifier score descending
    ///   - Walk down tracking cumulative targets (T) and decoys (D)
    ///   - FDR at each rank = D / max(T, 1)
    ///   - Q-value = monotonized minimum FDR from bottom up
    ///   - Populates DiaFdrInfo on each DiaSearchResult
    /// 
    /// Key design decision for iteration 1:
    ///   Uses DECOYS as negatives (not bottom-40% targets). This is critical because:
    ///   - Pseudo-labels from target quantiles only separate good-from-bad targets
    ///   - The classifier must learn what distinguishes targets FROM decoys
    ///   - Decoys may have pathological feature values (e.g., missing RT → extreme RtDeviation)
    ///     that the classifier must learn to penalize
    /// 
    /// Lives in MassSpectrometry/Dia/Scoring/ alongside other scoring classes.
    /// </summary>
    public static class DiaFdrEngine
    {
        // ════════════════════════════════════════════════════════════════
        //  Result Types
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Per-iteration diagnostics for monitoring convergence and target-decoy separation.
        /// </summary>
        public readonly struct IterationDiagnostics
        {
            public readonly int Iteration;
            public readonly int PositiveTrainingCount;
            public readonly int NegativeTrainingCount;
            public readonly float TargetScoreMean;
            public readonly float TargetScoreMedian;
            public readonly float TargetScoreQ25;
            public readonly float TargetScoreQ75;
            public readonly float DecoyScoreMean;
            public readonly float DecoyScoreMedian;
            public readonly float DecoyScoreQ25;
            public readonly float DecoyScoreQ75;
            public readonly float Separation;
            public readonly double Auc;
            public readonly int IdentificationsAt1Pct;
            public readonly float WeightChange;
            public readonly float[] Weights;
            public readonly float Bias;

            public IterationDiagnostics(
                int iteration, int positiveCount, int negativeCount,
                float targetMean, float targetMedian, float targetQ25, float targetQ75,
                float decoyMean, float decoyMedian, float decoyQ25, float decoyQ75,
                float separation, double auc, int idsAt1Pct,
                float weightChange, float[] weights, float bias)
            {
                Iteration = iteration;
                PositiveTrainingCount = positiveCount;
                NegativeTrainingCount = negativeCount;
                TargetScoreMean = targetMean;
                TargetScoreMedian = targetMedian;
                TargetScoreQ25 = targetQ25;
                TargetScoreQ75 = targetQ75;
                DecoyScoreMean = decoyMean;
                DecoyScoreMedian = decoyMedian;
                DecoyScoreQ25 = decoyQ25;
                DecoyScoreQ75 = decoyQ75;
                Separation = separation;
                Auc = auc;
                IdentificationsAt1Pct = idsAt1Pct;
                WeightChange = weightChange;
                Weights = weights;
                Bias = bias;
            }
        }

        /// <summary>
        /// Result of the full iterative FDR estimation pipeline.
        /// </summary>
        public readonly struct FdrResult
        {
            public readonly DiaLinearDiscriminant Classifier;
            public readonly int IterationsCompleted;
            public readonly int IdentificationsAt1PctFdr;
            public readonly IterationDiagnostics[] Diagnostics;

            public FdrResult(DiaLinearDiscriminant classifier, int iterations,
                int idsAt1Pct, IterationDiagnostics[] diagnostics)
            {
                Classifier = classifier;
                IterationsCompleted = iterations;
                IdentificationsAt1PctFdr = idsAt1Pct;
                Diagnostics = diagnostics;
            }
        }

        // ════════════════════════════════════════════════════════════════
        //  Main Entry Point
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Runs the full iterative target-decoy FDR estimation pipeline.
        /// 
        /// Modifies ClassifierScore and FdrInfo on each DiaSearchResult in place.
        /// Features array must be parallel to results (same indices).
        /// </summary>
        public static FdrResult RunIterativeFdr(
            List<DiaSearchResult> results,
            DiaFeatureVector[] features,
            int maxIterations = 5,
            float convergenceThreshold = 0.01f,
            float idCountConvergenceThreshold = 0.01f,
            float l2Lambda = 5e-3f,
            float learningRate = 0.05f,
            int maxEpochs = 300)
        {
            if (results == null || results.Count == 0)
                throw new ArgumentException("Results list cannot be null or empty.");
            if (features == null || features.Length != results.Count)
                throw new ArgumentException("Features array must be parallel to results.");

            var diagnosticsList = new List<IterationDiagnostics>();
            float[] previousWeights = null;
            int previousIdCount = 0;
            DiaLinearDiscriminant classifier = null;

            // ── Separate target and decoy indices ───────────────────────
            var targetIndices = new List<int>(results.Count);
            var decoyIndices = new List<int>(results.Count);
            for (int i = 0; i < results.Count; i++)
            {
                if (results[i].IsDecoy)
                    decoyIndices.Add(i);
                else
                    targetIndices.Add(i);
            }

            if (decoyIndices.Count == 0)
                throw new InvalidOperationException("No decoy results found. Cannot compute FDR.");
            if (targetIndices.Count == 0)
                throw new InvalidOperationException("No target results found.");

            for (int iteration = 1; iteration <= maxIterations; iteration++)
            {
                DiaFeatureVector[] positives;
                DiaFeatureVector[] negatives;

                if (iteration == 1)
                {
                    // ── Iteration 1: top 40% targets as POSITIVES, decoys as NEGATIVES ──
                    // Using decoys as negatives (not bottom-40% targets) is critical:
                    // the classifier must learn what distinguishes targets from decoys,
                    // not just good targets from bad targets.
                    positives = BuildPseudoLabelPositives(features, targetIndices);
                    negatives = CollectDecoyNegatives(features, decoyIndices, positives.Length);
                }
                else
                {
                    // ── Iteration 2+: confident targets vs decoys ───────────
                    positives = CollectConfidentTargets(results, features, targetIndices, qThreshold: 0.01f);
                    negatives = CollectDecoyNegatives(features, decoyIndices, positives.Length);
                }

                if (positives.Length < 50 || negatives.Length < 50)
                    break;

                // ── Train classifier ────────────────────────────────────
                classifier = DiaLinearDiscriminant.TrainLogisticRegression(
                    positives.AsSpan(), negatives.AsSpan(),
                    learningRate: learningRate,
                    l2Lambda: l2Lambda,
                    maxEpochs: maxEpochs,
                    batchSize: 256,
                    useInteraction: false);

                // ── Score all results ───────────────────────────────────
                for (int i = 0; i < results.Count; i++)
                    results[i].ClassifierScore = classifier.Score(in features[i]);

                // ── Compute q-values ────────────────────────────────────
                int idsAt1Pct = ComputeQValues(results);

                // ── Diagnostics ─────────────────────────────────────────
                float weightChange = previousWeights != null
                    ? classifier.WeightChangeL2(previousWeights)
                    : float.PositiveInfinity;

                var diag = BuildDiagnostics(
                    iteration, results, targetIndices, decoyIndices,
                    positives.Length, negatives.Length,
                    idsAt1Pct, weightChange, classifier);
                diagnosticsList.Add(diag);

                // ── Convergence check (skip iteration 1) ────────────────
                if (iteration > 1)
                {
                    bool converged = false;

                    if (weightChange < convergenceThreshold)
                        converged = true;

                    if (previousIdCount > 0)
                    {
                        float idChangeRatio = MathF.Abs(idsAt1Pct - previousIdCount)
                                              / (float)previousIdCount;
                        if (idChangeRatio < idCountConvergenceThreshold)
                            converged = true;
                    }

                    if (converged)
                    {
                        previousWeights = (float[])classifier.Weights.Clone();
                        previousIdCount = idsAt1Pct;
                        break;
                    }
                }

                previousWeights = (float[])classifier.Weights.Clone();
                previousIdCount = idsAt1Pct;
            }

            // ── Final pass ──────────────────────────────────────────────
            int finalIds = 0;
            if (classifier != null)
            {
                for (int i = 0; i < results.Count; i++)
                    results[i].ClassifierScore = classifier.Score(in features[i]);
                finalIds = ComputeQValues(results);
            }

            return new FdrResult(
                classifier,
                diagnosticsList.Count,
                finalIds,
                diagnosticsList.ToArray());
        }

        // ════════════════════════════════════════════════════════════════
        //  Q-Value Computation
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Computes q-values for all results using the target-decoy approach.
        /// 
        /// Algorithm:
        ///   1. Sort by ClassifierScore descending
        ///   2. Walk down tracking cumulative T and D
        ///   3. FDR = D / max(T, 1)
        ///   4. Monotonize from bottom up (q-value = running min)
        ///   5. Populate DiaFdrInfo on each result
        /// 
        /// Performance: O(n log n) for sort + O(n) for forward/backward passes.
        /// Returns count of target identifications at q ≤ 0.01.
        /// </summary>
        public static int ComputeQValues(List<DiaSearchResult> results)
        {
            int n = results.Count;
            if (n == 0) return 0;

            int[] sortedIdx = ArrayPool<int>.Shared.Rent(n);
            double[] rawFdr = ArrayPool<double>.Shared.Rent(n);
            int[] cumTargets = ArrayPool<int>.Shared.Rent(n);
            int[] cumDecoys = ArrayPool<int>.Shared.Rent(n);
            try
            {
                for (int i = 0; i < n; i++)
                    sortedIdx[i] = i;

                // Sort descending by ClassifierScore
                Array.Sort(sortedIdx, 0, n, new DescendingScoreComparer(results));

                // Forward pass: compute cumulative T, D, and raw FDR in single O(n) sweep
                int t = 0, d = 0;
                for (int rank = 0; rank < n; rank++)
                {
                    int idx = sortedIdx[rank];
                    if (results[idx].IsDecoy)
                        d++;
                    else
                        t++;

                    cumTargets[rank] = t;
                    cumDecoys[rank] = d;
                    rawFdr[rank] = t > 0 ? (double)d / t : 1.0;
                }

                // Backward pass: monotonize q-values (running minimum from bottom)
                double runningMin = 1.0;
                for (int rank = n - 1; rank >= 0; rank--)
                {
                    if (rawFdr[rank] < runningMin)
                        runningMin = rawFdr[rank];
                    rawFdr[rank] = runningMin;
                }

                // Assign to results — O(n), no re-counting
                int idsAt1Pct = 0;
                for (int rank = 0; rank < n; rank++)
                {
                    int idx = sortedIdx[rank];
                    var result = results[idx];

                    if (result.FdrInfo == null)
                        result.FdrInfo = new DiaFdrInfo();

                    result.FdrInfo.CumulativeTarget = cumTargets[rank];
                    result.FdrInfo.CumulativeDecoy = cumDecoys[rank];
                    result.FdrInfo.QValue = Math.Min(rawFdr[rank], 1.0);

                    if (!result.IsDecoy && result.FdrInfo.QValue <= 0.01)
                        idsAt1Pct++;
                }

                return idsAt1Pct;
            }
            finally
            {
                ArrayPool<int>.Shared.Return(sortedIdx);
                ArrayPool<double>.Shared.Return(rawFdr);
                ArrayPool<int>.Shared.Return(cumTargets);
                ArrayPool<int>.Shared.Return(cumDecoys);
            }
        }

        /// <summary>Comparer for descending sort by ClassifierScore.</summary>
        private sealed class DescendingScoreComparer : IComparer<int>
        {
            private readonly List<DiaSearchResult> _results;
            public DescendingScoreComparer(List<DiaSearchResult> results) => _results = results;
            public int Compare(int a, int b) =>
                _results[b].ClassifierScore.CompareTo(_results[a].ClassifierScore);
        }

        // ════════════════════════════════════════════════════════════════
        //  Training Data Construction
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Top 40% of targets by max(ApexScore, TemporalScore) as positive training set.
        /// </summary>
        private static DiaFeatureVector[] BuildPseudoLabelPositives(
            DiaFeatureVector[] features, List<int> targetIndices)
        {
            int nTargets = targetIndices.Count;
            float[] scores = ArrayPool<float>.Shared.Rent(nTargets);
            int[] order = ArrayPool<int>.Shared.Rent(nTargets);
            try
            {
                for (int i = 0; i < nTargets; i++)
                {
                    int idx = targetIndices[i];
                    scores[i] = MathF.Max(features[idx].ApexScore, features[idx].TemporalScore);
                    order[i] = i;
                }
                Array.Sort(scores, order, 0, nTargets);

                int cutoff = (int)(nTargets * 0.6);
                int posCount = nTargets - cutoff;
                var positives = new DiaFeatureVector[posCount];
                for (int i = 0; i < posCount; i++)
                    positives[i] = features[targetIndices[order[cutoff + i]]];
                return positives;
            }
            finally
            {
                ArrayPool<float>.Shared.Return(scores);
                ArrayPool<int>.Shared.Return(order);
            }
        }

        /// <summary>
        /// Collects confident target identifications (q &lt; threshold) as positives.
        /// </summary>
        private static DiaFeatureVector[] CollectConfidentTargets(
            List<DiaSearchResult> results, DiaFeatureVector[] features,
            List<int> targetIndices, float qThreshold)
        {
            var positives = new List<DiaFeatureVector>(targetIndices.Count / 2);
            for (int i = 0; i < targetIndices.Count; i++)
            {
                int idx = targetIndices[i];
                var fdr = results[idx].FdrInfo;
                if (fdr != null && fdr.QValue <= qThreshold)
                    positives.Add(features[idx]);
            }
            return positives.ToArray();
        }

        /// <summary>
        /// Collects decoy feature vectors as negatives. Subsamples to balance if needed.
        /// </summary>
        private static DiaFeatureVector[] CollectDecoyNegatives(
            DiaFeatureVector[] features, List<int> decoyIndices, int positiveCount)
        {
            if (decoyIndices.Count <= positiveCount * 2)
            {
                var negatives = new DiaFeatureVector[decoyIndices.Count];
                for (int i = 0; i < decoyIndices.Count; i++)
                    negatives[i] = features[decoyIndices[i]];
                return negatives;
            }

            // Subsample decoys to match positive count
            var rng = new Random(42);
            int sampleSize = Math.Max(positiveCount, 1000); // ensure minimum diversity
            sampleSize = Math.Min(sampleSize, decoyIndices.Count);
            var sampled = new DiaFeatureVector[sampleSize];
            var shuffled = new int[decoyIndices.Count];
            for (int i = 0; i < decoyIndices.Count; i++)
                shuffled[i] = decoyIndices[i];

            for (int i = 0; i < sampleSize; i++)
            {
                int j = i + rng.Next(shuffled.Length - i);
                (shuffled[i], shuffled[j]) = (shuffled[j], shuffled[i]);
                sampled[i] = features[shuffled[i]];
            }
            return sampled;
        }

        // ════════════════════════════════════════════════════════════════
        //  Diagnostics
        // ════════════════════════════════════════════════════════════════

        private static IterationDiagnostics BuildDiagnostics(
            int iteration,
            List<DiaSearchResult> results,
            List<int> targetIndices,
            List<int> decoyIndices,
            int positiveTrainCount,
            int negativeTrainCount,
            int idsAt1Pct,
            float weightChange,
            DiaLinearDiscriminant classifier)
        {
            float[] targetScores = new float[targetIndices.Count];
            for (int i = 0; i < targetIndices.Count; i++)
                targetScores[i] = results[targetIndices[i]].ClassifierScore;

            float[] decoyScores = new float[decoyIndices.Count];
            for (int i = 0; i < decoyIndices.Count; i++)
                decoyScores[i] = results[decoyIndices[i]].ClassifierScore;

            Array.Sort(targetScores);
            Array.Sort(decoyScores);

            var tStats = QuartileStats(targetScores);
            var dStats = QuartileStats(decoyScores);
            float separation = tStats.Mean - dStats.Mean;

            double auc = ComputeTargetDecoyAuc(results, targetIndices, decoyIndices);
            float[] weightsCopy = (float[])classifier.Weights.Clone();

            return new IterationDiagnostics(
                iteration, positiveTrainCount, negativeTrainCount,
                tStats.Mean, tStats.Median, tStats.Q25, tStats.Q75,
                dStats.Mean, dStats.Median, dStats.Q25, dStats.Q75,
                separation, auc, idsAt1Pct,
                weightChange, weightsCopy, classifier.Bias);
        }

        private readonly struct Stats
        {
            public readonly float Mean, Median, Q25, Q75;
            public Stats(float mean, float median, float q25, float q75)
            { Mean = mean; Median = median; Q25 = q25; Q75 = q75; }
        }

        private static Stats QuartileStats(float[] sorted)
        {
            if (sorted.Length == 0) return new Stats(0, 0, 0, 0);
            int n = sorted.Length;
            float mean = 0;
            for (int i = 0; i < n; i++) mean += sorted[i];
            mean /= n;
            return new Stats(mean, sorted[n / 2], sorted[n / 4], sorted[3 * n / 4]);
        }

        private static double ComputeTargetDecoyAuc(
            List<DiaSearchResult> results,
            List<int> targetIndices,
            List<int> decoyIndices)
        {
            int totalPos = targetIndices.Count;
            int totalNeg = decoyIndices.Count;
            if (totalPos == 0 || totalNeg == 0) return 0.5;

            var scored = new (float Score, bool IsTarget)[totalPos + totalNeg];
            for (int i = 0; i < totalPos; i++)
                scored[i] = (results[targetIndices[i]].ClassifierScore, true);
            for (int i = 0; i < totalNeg; i++)
                scored[totalPos + i] = (results[decoyIndices[i]].ClassifierScore, false);

            Array.Sort(scored, (a, b) => b.Score.CompareTo(a.Score));

            double auc = 0, prevFpr = 0, prevTpr = 0;
            int tp = 0, fp = 0;
            for (int i = 0; i < scored.Length; i++)
            {
                if (scored[i].IsTarget) tp++;
                else fp++;
                double tpr = (double)tp / totalPos;
                double fpr = (double)fp / totalNeg;
                auc += (fpr - prevFpr) * (tpr + prevTpr) / 2.0;
                prevFpr = fpr;
                prevTpr = tpr;
            }
            auc += (1.0 - prevFpr) * (1.0 + prevTpr) / 2.0;
            return auc;
        }

        // ════════════════════════════════════════════════════════════════
        //  Console Reporting
        // ════════════════════════════════════════════════════════════════

        public static void PrintDiagnostics(IterationDiagnostics d)
        {
            Console.WriteLine($"  ─── Iteration {d.Iteration} ───");
            Console.WriteLine($"    Training: {d.PositiveTrainingCount:N0} positives, {d.NegativeTrainingCount:N0} negatives");
            Console.WriteLine($"    Target scores:  mean={d.TargetScoreMean:F4}  median={d.TargetScoreMedian:F4}  Q25={d.TargetScoreQ25:F4}  Q75={d.TargetScoreQ75:F4}");
            Console.WriteLine($"    Decoy scores:   mean={d.DecoyScoreMean:F4}  median={d.DecoyScoreMedian:F4}  Q25={d.DecoyScoreQ25:F4}  Q75={d.DecoyScoreQ75:F4}");
            Console.WriteLine($"    Separation:     {d.Separation:F4} (mean_target - mean_decoy)");
            Console.WriteLine($"    AUC:            {d.Auc:F4}");
            Console.WriteLine($"    IDs at 1% FDR:  {d.IdentificationsAt1Pct:N0}");
            Console.WriteLine($"    Weight change:  {(float.IsPositiveInfinity(d.WeightChange) ? "∞" : d.WeightChange.ToString("F6"))}");
            Console.WriteLine($"    Bias:           {d.Bias:F6}");

            int nWeights = Math.Min(d.Weights.Length, DiaFeatureVector.ClassifierFeatureCount);
            var indexed = new (string Name, float Weight)[nWeights];
            for (int i = 0; i < nWeights; i++)
                indexed[i] = (DiaFeatureVector.FeatureNames[i], d.Weights[i]);
            Array.Sort(indexed, (a, b) => MathF.Abs(b.Weight).CompareTo(MathF.Abs(a.Weight)));
            Console.WriteLine("    Weights:");
            for (int i = 0; i < nWeights; i++)
            {
                string sign = indexed[i].Weight >= 0 ? "+" : "";
                Console.WriteLine($"      {sign}{indexed[i].Weight:F4}  {indexed[i].Name}");
            }
        }
    }
}