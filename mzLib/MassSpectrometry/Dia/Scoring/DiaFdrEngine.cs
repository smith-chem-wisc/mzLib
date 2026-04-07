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
    /// Implements the iterative semi-supervised learning approach:
    ///   Iteration 1: Train on pseudo-labels — top 40% targets (positive) vs ALL decoys (negative).
    ///   Iteration 2+: Retrain using q &lt; 0.01 targets (positive) vs subsampled decoys (negative).
    ///   Repeat until convergence (weight change &lt; 1%, or ID count stable, or max iterations).
    /// 
    /// Phase 14 changes:
    ///   - IDiaClassifier abstraction: classifiers are swappable at runtime
    ///   - DiaGradientBoostedClassifier (GBT) as the default non-linear classifier
    ///   - DiaLinearDiscriminant (LDA) preserved via DiaLdaClassifierAdapter
    ///   - RunIterativeFdr() gains a DiaClassifierType parameter (default: GBT)
    ///   - Diagnostics extended with classifier type info
    ///   - GBT feature importances reported alongside LDA weights
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
            public readonly DiaClassifierType ClassifierType;

            // LDA-specific (null for GBT)
            public readonly float[] Weights;
            public readonly float Bias;

            // GBT-specific (null for LDA)
            public readonly float[] FeatureImportances;

            public IterationDiagnostics(
                int iteration, int positiveCount, int negativeCount,
                float targetMean, float targetMedian, float targetQ25, float targetQ75,
                float decoyMean, float decoyMedian, float decoyQ25, float decoyQ75,
                float separation, double auc, int idsAt1Pct,
                float weightChange, DiaClassifierType classifierType,
                float[] weights, float bias, float[] featureImportances)
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
                ClassifierType = classifierType;
                Weights = weights;
                Bias = bias;
                FeatureImportances = featureImportances;
            }
        }

        /// <summary>
        /// Result of the full iterative FDR estimation pipeline.
        /// </summary>
        public readonly struct FdrResult
        {
            /// <summary>The trained classifier (LDA adapter or GBT adapter)</summary>
            public readonly IDiaClassifier Classifier;

            /// <summary>If LDA was used, the underlying DiaLinearDiscriminant. Null for GBT.</summary>
            public readonly DiaLinearDiscriminant LdaClassifier;

            public readonly DiaClassifierType ClassifierType;
            public readonly int IterationsCompleted;
            public readonly int IdentificationsAt1PctFdr;
            public readonly IterationDiagnostics[] Diagnostics;

            public FdrResult(IDiaClassifier classifier, DiaLinearDiscriminant ldaClassifier,
                DiaClassifierType classifierType, int iterations,
                int idsAt1Pct, IterationDiagnostics[] diagnostics)
            {
                Classifier = classifier;
                LdaClassifier = ldaClassifier;
                ClassifierType = classifierType;
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
        /// <param name="results">All search results (targets + decoys). Modified in-place.</param>
        /// <param name="features">Feature vectors parallel to results.</param>
        /// <param name="classifierType">LDA (Phase 10) or GBT (Phase 14, default).</param>
        /// <param name="maxIterations">Maximum semi-supervised iterations.</param>
        /// <param name="convergenceThreshold">LDA weight change threshold for early stop.</param>
        /// <param name="idCountConvergenceThreshold">ID count change ratio threshold for early stop.</param>
        /// <param name="l2Lambda">L2 regularization for LDA. Ignored for GBT.</param>
        /// <param name="learningRate">Learning rate for LDA. Ignored for GBT (uses its own).</param>
        /// <param name="maxEpochs">Max SGD epochs for LDA. Ignored for GBT.</param>
        public static FdrResult RunIterativeFdr(
            List<DiaSearchResult> results,
            DiaFeatureVector[] features,
            DiaClassifierType classifierType = DiaClassifierType.GradientBoostedTree,
            int maxIterations = 10,
            float convergenceThreshold = 0.01f,
            float idCountConvergenceThreshold = 0.005f,
            float l2Lambda = 5e-3f,
            float learningRate = 0.05f,
            int maxEpochs = 300,
            bool rtWindowRefinement = false,
            float rtSigmaMultiplier = 3.0f,
            Action<List<DiaSearchResult>, DiaFeatureVector[]> afterIteration1 = null)
        {
            if (results == null || results.Count == 0)
                throw new ArgumentException("Results list cannot be null or empty.");
            if (features == null || features.Length != results.Count)
                throw new ArgumentException("Features array must be parallel to results.");

            var diagnosticsList = new List<IterationDiagnostics>();
            int previousIdCount = 0;
            IDiaClassifier classifier = null;

            // ── Best-iteration tracking ─────────────────────────────────
            // GBT can oscillate: iteration 2 finds 14K IDs, iteration 3 drops to 12K.
            // We keep the classifier that produced the most IDs and restore it at the end.
            IDiaClassifier bestClassifier = null;
            int bestIds = 0;
            int bestIteration = 0;

            // For LDA weight-change convergence tracking
            float[] previousLdaWeights = null;

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

            // Track consecutive declines to detect oscillation
            int consecutiveDeclines = 0;

            // ── RT window refinement state ──────────────────────────────
            // active[i] = false means row i is excluded from training and FDR
            // for this iteration (its RT deviation is outside the current window).
            // We never permanently discard rows — the mask is recomputed each iteration.
            bool[] active = new bool[results.Count];
            for (int i = 0; i < active.Length; i++) active[i] = true;
            float rtSigma = float.MaxValue;

            for (int iteration = 1; iteration <= maxIterations; iteration++)
            {
                DiaFeatureVector[] positives;
                DiaFeatureVector[] negatives;

                if (iteration == 1)
                {
                    // ── Iteration 1: top 40% targets as POSITIVES, decoys as NEGATIVES ──
                    // All rows are active on the first pass — we have no estimate yet.
                    positives = BuildPseudoLabelPositives(features, targetIndices, active);
                    negatives = CollectDecoyNegatives(features, decoyIndices, positives.Length, active);
                }
                else
                {
                    // ── Iteration 2+: confident targets vs decoys ───────────
                    positives = CollectConfidentTargets(results, features, targetIndices, qThreshold: 0.01f, active);
                    negatives = CollectDecoyNegatives(features, decoyIndices, positives.Length, active);
                }

                if (positives.Length < 50 || negatives.Length < 50)
                    break;

                // ── Train classifier ────────────────────────────────────
                classifier = CreateClassifier(classifierType, learningRate, l2Lambda, maxEpochs);
                classifier.Train(positives.AsSpan(), negatives.AsSpan());

                // ── Score all results ───────────────────────────────────
                for (int i = 0; i < results.Count; i++)
                    results[i].ClassifierScore = classifier.Score(in features[i]);

                // ── Compute q-values (active rows only) ─────────────────
                int idsAt1Pct = ComputeQValues(results, active);

                // ── RT window refinement (local bias correction) ───────────
                // Build a local bias curve from confident targets: for each result,
                // the "center" is the median deviation of confident neighbors within
                // ±localWindowMinutes of its observed apex RT. This corrects the
                // systematic tilt in the LOWESS calibration without re-fitting it.
                // The window filter is: |dev - localBias| ≤ rtSigmaMultiplier * sigma
                // where sigma measures spread around the local center (not global).
                if (rtWindowRefinement && iteration >= 1)
                {
                    float[] localBias;
                    float newSigma;
                    const float localWindowMinutes = 2.0f;  // RT neighborhood for bias estimate

                    ComputeLocalRtBias(results, targetIndices, active,
                        localWindowMinutes, out localBias, out newSigma);

                    if (!float.IsNaN(newSigma) && newSigma > 0f)
                    {
                        float hw = rtSigmaMultiplier * newSigma;
                        int activeBefore = 0, activeAfter = 0;
                        for (int i = 0; i < active.Length; i++) if (active[i]) activeBefore++;

                        for (int i = 0; i < results.Count; i++)
                        {
                            float dev = results[i].RtDeviationMinutes;
                            float bias = localBias[i];
                            // If we have no local bias estimate (sparse region), fall back
                            // to the raw deviation with a wider tolerance.
                            active[i] = !float.IsNaN(dev) && (
                                float.IsNaN(bias)
                                    ? MathF.Abs(dev) <= hw * 1.5f   // fallback: wider window
                                    : MathF.Abs(dev - bias) <= hw);
                        }

                        for (int i = 0; i < active.Length; i++) if (active[i]) activeAfter++;

                        int activeT = 0, activeD = 0;
                        for (int i = 0; i < results.Count; i++)
                        {
                            if (!active[i]) continue;
                            if (results[i].IsDecoy) activeD++; else activeT++;
                        }
                        Console.WriteLine($"  [RTrefine] iter={iteration}  " +
                            $"localWindow=±{localWindowMinutes:F1} min  " +
                            $"σ={newSigma:F3} min  " +
                            $"hw={hw:F3} min  " +
                            $"active: {activeBefore:N0}→{activeAfter:N0}  " +
                            $"T/D={activeT:N0}/{activeD:N0}≈{(activeD > 0 ? (float)activeT / activeD : 0f):F2}");

                        rtSigma = newSigma;
                    }
                }

                // ── Track best iteration ────────────────────────────────
                if (idsAt1Pct > bestIds)
                {
                    bestIds = idsAt1Pct;
                    bestClassifier = classifier;
                    bestIteration = iteration;
                    consecutiveDeclines = 0;
                }
                else
                {
                    consecutiveDeclines++;
                }

                // ── Weight change (LDA only; GBT uses ID count convergence) ──
                float weightChange = float.PositiveInfinity;
                if (classifierType == DiaClassifierType.LinearDiscriminant &&
                    classifier is DiaLdaClassifierAdapter ldaAdapter && ldaAdapter.Lda != null)
                {
                    if (previousLdaWeights != null)
                        weightChange = ldaAdapter.Lda.WeightChangeL2(previousLdaWeights);
                    previousLdaWeights = (float[])ldaAdapter.Lda.Weights.Clone();
                }

                // ── Diagnostics ─────────────────────────────────────────
                var diag = BuildDiagnostics(
                    iteration, results, targetIndices, decoyIndices,
                    positives.Length, negatives.Length,
                    idsAt1Pct, weightChange, classifierType, classifier);
                diagnosticsList.Add(diag);

                // ── Post-iteration-1 injection hook ─────────────────────
                // Fires once after iteration 1 has written ClassifierScore
                // onto all results. Allows caller to inject features that
                // depend on relative classifier scores across results
                // (e.g. ChargeStateRtConsensus). The updated feature vectors
                // are then used when iteration 2 collects training positives.
                if (iteration == 1 && afterIteration1 != null)
                    afterIteration1(results, features);

                // ── Convergence check (skip iteration 1) ────────────────
                if (iteration > 1)
                {
                    bool converged = false;

                    // LDA: weight change convergence
                    if (classifierType == DiaClassifierType.LinearDiscriminant &&
                        weightChange < convergenceThreshold)
                        converged = true;

                    // Both: ID count convergence (stable)
                    if (previousIdCount > 0)
                    {
                        float idChangeRatio = MathF.Abs(idsAt1Pct - previousIdCount)
                                              / (float)previousIdCount;
                        if (idChangeRatio < idCountConvergenceThreshold)
                            converged = true;
                    }

                    // GBT: stop if declining for 2 consecutive iterations after peak
                    // This catches the oscillation pattern: 12K → 14K → 12K → 12K
                    if (classifierType == DiaClassifierType.GradientBoostedTree &&
                        consecutiveDeclines >= 2 && bestIds > 0)
                        converged = true;

                    if (converged)
                    {
                        previousIdCount = idsAt1Pct;
                        break;
                    }
                }

                previousIdCount = idsAt1Pct;
            }

            // ── Final pass: restore best classifier ─────────────────────
            // For GBT, the best iteration often isn't the last one due to oscillation.
            // Re-score with the best classifier and recompute q-values.
            int finalIds = 0;
            IDiaClassifier finalClassifier = bestClassifier ?? classifier;
            if (finalClassifier != null)
            {
                if (finalClassifier != classifier)
                {
                    // Best wasn't last — re-score everything with the best model
                    for (int i = 0; i < results.Count; i++)
                        results[i].ClassifierScore = finalClassifier.Score(in features[i]);
                }
                finalIds = ComputeQValues(results, active);

                if (bestIteration > 0 && bestClassifier != classifier)
                {
                    Console.WriteLine($"  [Best iteration: {bestIteration} with {bestIds:N0} IDs — restored]");
                }
            }

            // Extract LDA classifier if applicable
            DiaLinearDiscriminant ldaOut = null;
            if (finalClassifier is DiaLdaClassifierAdapter finalLda)
                ldaOut = finalLda.Lda;

            return new FdrResult(
                finalClassifier,
                ldaOut,
                classifierType,
                diagnosticsList.Count,
                finalIds,
                diagnosticsList.ToArray());
        }

        // ════════════════════════════════════════════════════════════════
        //  Classifier Factory
        // ════════════════════════════════════════════════════════════════

        private static IDiaClassifier CreateClassifier(
            DiaClassifierType type,
            float learningRate, float l2Lambda, int maxEpochs)
        {
            return type switch
            {
                DiaClassifierType.GradientBoostedTree => new DiaGbtClassifierAdapter(
                    new DiaGradientBoostedClassifier(
                        featureCount: DiaFeatureVector.ClassifierFeatureCount,
                        numTrees: 200,
                        maxDepth: 5,
                        learningRate: 0.05f,
                        minSamplesLeaf: 15,
                        l2Regularization: 0.5f,
                        numBins: 128,
                        subsampleFraction: 0.85f)),

                DiaClassifierType.LinearDiscriminant => new DiaLdaClassifierAdapter(
                    learningRate: learningRate,
                    l2Lambda: l2Lambda,
                    maxEpochs: maxEpochs,
                    batchSize: 256),

                DiaClassifierType.NeuralNetwork => new DiaNeuralNetClassifierAdapter(),

                _ => throw new ArgumentException($"Unknown classifier type: {type}")
            };
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
        /// Returns count of target identifications at q ≤ 0.01.
        /// </summary>
        public static int ComputeQValues(List<DiaSearchResult> results, bool[] active = null)
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

                Array.Sort(sortedIdx, 0, n, new DescendingScoreComparer(results));

                int t = 0, d = 0;
                for (int rank = 0; rank < n; rank++)
                {
                    int idx = sortedIdx[rank];
                    // Inactive rows are scored but assigned maximum FDR (treated as noise).
                    // They do not contribute to cumulative T/D counts.
                    bool isActive = active == null || active[idx];
                    if (isActive)
                    {
                        if (results[idx].IsDecoy) d++;
                        else t++;
                    }

                    cumTargets[rank] = t;
                    cumDecoys[rank] = d;
                    rawFdr[rank] = t > 0 ? (double)d / t : 1.0;
                }

                double runningMin = 1.0;
                for (int rank = n - 1; rank >= 0; rank--)
                {
                    if (rawFdr[rank] < runningMin)
                        runningMin = rawFdr[rank];
                    rawFdr[rank] = runningMin;
                }

                int idsAt1Pct = 0;
                for (int rank = 0; rank < n; rank++)
                {
                    int idx = sortedIdx[rank];
                    var result = results[idx];

                    if (result.FdrInfo == null)
                        result.FdrInfo = new DiaFdrInfo();

                    result.FdrInfo.CumulativeTarget = cumTargets[rank];
                    result.FdrInfo.CumulativeDecoy = cumDecoys[rank];

                    // Inactive rows get q-value = 1.0 regardless of score.
                    bool isActive = active == null || active[idx];
                    result.FdrInfo.QValue = isActive ? Math.Min(rawFdr[rank], 1.0) : 1.0;

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
            DiaFeatureVector[] features, List<int> targetIndices, bool[] active)
        {
            int nTargets = targetIndices.Count;
            float[] scores = ArrayPool<float>.Shared.Rent(nTargets);
            int[] order = ArrayPool<int>.Shared.Rent(nTargets);
            try
            {
                for (int i = 0; i < nTargets; i++)
                {
                    int idx = targetIndices[i];
                    scores[i] = active[idx]
                        ? MathF.Max(features[idx].ApexScore, features[idx].TemporalScore)
                        : float.MinValue;  // inactive: sort to bottom
                    order[i] = i;
                }
                Array.Sort(scores, order, 0, nTargets);

                int cutoff = (int)(nTargets * 0.6);
                int posCount = nTargets - cutoff;
                var positives = new List<DiaFeatureVector>(posCount);
                for (int i = cutoff; i < nTargets; i++)
                {
                    int idx = targetIndices[order[i]];
                    if (active[idx])
                        positives.Add(features[idx]);
                }
                return positives.ToArray();
            }
            finally
            {
                ArrayPool<float>.Shared.Return(scores);
                ArrayPool<int>.Shared.Return(order);
            }
        }

        private static DiaFeatureVector[] CollectConfidentTargets(
            List<DiaSearchResult> results, DiaFeatureVector[] features,
            List<int> targetIndices, float qThreshold, bool[] active)
        {
            var positives = new List<DiaFeatureVector>(targetIndices.Count / 2);
            for (int i = 0; i < targetIndices.Count; i++)
            {
                int idx = targetIndices[i];
                if (!active[idx]) continue;
                var fdr = results[idx].FdrInfo;
                if (fdr != null && fdr.QValue <= qThreshold)
                    positives.Add(features[idx]);
            }
            return positives.ToArray();
        }

        /// <summary>
        /// Collects decoy feature vectors as negatives. Subsamples to balance if needed.
        /// Only includes active rows.
        /// </summary>
        private static DiaFeatureVector[] CollectDecoyNegatives(
            DiaFeatureVector[] features, List<int> decoyIndices, int positiveCount, bool[] active)
        {
            // Count active decoys
            int activeDecoyCount = 0;
            for (int i = 0; i < decoyIndices.Count; i++)
                if (active[decoyIndices[i]]) activeDecoyCount++;

            if (activeDecoyCount <= positiveCount * 2)
            {
                var negatives = new DiaFeatureVector[activeDecoyCount];
                int k = 0;
                for (int i = 0; i < decoyIndices.Count; i++)
                {
                    int idx = decoyIndices[i];
                    if (active[idx]) negatives[k++] = features[idx];
                }
                return negatives;
            }

            // Subsample active decoys
            var rng = new Random(42);
            int sampleSize = Math.Max(positiveCount, 1000);
            sampleSize = Math.Min(sampleSize, activeDecoyCount);

            // Collect active decoy indices then shuffle-sample
            var activeDecoys = new int[activeDecoyCount];
            int j = 0;
            for (int i = 0; i < decoyIndices.Count; i++)
            {
                int idx = decoyIndices[i];
                if (active[idx]) activeDecoys[j++] = idx;
            }

            var sampled = new DiaFeatureVector[sampleSize];
            for (int i = 0; i < sampleSize; i++)
            {
                int swap = i + rng.Next(activeDecoys.Length - i);
                (activeDecoys[i], activeDecoys[swap]) = (activeDecoys[swap], activeDecoys[i]);
                sampled[i] = features[activeDecoys[i]];
            }
            return sampled;
        }

        // ════════════════════════════════════════════════════════════════
        //  RT Window Refinement Helpers
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Builds a local bias correction curve from confident targets.
        /// 
        /// For each result, the "local bias" is the median RtDeviationMinutes of
        /// confident targets (q ≤ 0.01) whose ObservedApexRt falls within
        /// ±localWindowMinutes. This corrects the systematic tilt in the LOWESS
        /// calibration curve without requiring a re-fit.
        /// 
        /// Returns:
        ///   localBias[i]  — the local calibration offset for result i (NaN if no neighbors)
        ///   sigma         — global std of (dev - localBias) across all confident targets
        ///                   (measures spread around the local center, not the global center)
        /// </summary>
        private static void ComputeLocalRtBias(
            List<DiaSearchResult> results,
            List<int> targetIndices,
            bool[] active,
            float localWindowMinutes,
            out float[] localBias,
            out float sigma)
        {
            localBias = new float[results.Count];
            for (int i = 0; i < localBias.Length; i++) localBias[i] = float.NaN;

            // Collect confident targets sorted by ObservedApexRt
            var anchors = new List<(float apexRt, float dev)>(targetIndices.Count / 4);
            for (int i = 0; i < targetIndices.Count; i++)
            {
                int idx = targetIndices[i];
                if (!active[idx]) continue;
                var r = results[idx];
                if (r.FdrInfo == null || r.FdrInfo.QValue > 0.01f) continue;
                float dev = r.RtDeviationMinutes;
                float rt = r.ObservedApexRt;
                if (float.IsNaN(dev) || float.IsNaN(rt)) continue;
                anchors.Add((rt, dev));
            }

            if (anchors.Count < 10) { sigma = float.NaN; return; }

            anchors.Sort((a, b) => a.apexRt.CompareTo(b.apexRt));

            float[] anchorRts = new float[anchors.Count];
            float[] anchorDevs = new float[anchors.Count];
            for (int i = 0; i < anchors.Count; i++)
            {
                anchorRts[i] = anchors[i].apexRt;
                anchorDevs[i] = anchors[i].dev;
            }

            // For each result, binary-search into anchorRts to find neighbors
            // within ±localWindowMinutes, then take their median deviation.
            var neighborBuf = new List<float>(64);
            for (int i = 0; i < results.Count; i++)
            {
                float rt = results[i].ObservedApexRt;
                if (float.IsNaN(rt)) continue;

                float lo = rt - localWindowMinutes;
                float hi = rt + localWindowMinutes;

                // Binary search for left boundary
                int left = LowerBound(anchorRts, lo);
                int right = UpperBound(anchorRts, hi);

                if (right - left < 3) continue;  // too few neighbors — leave NaN

                neighborBuf.Clear();
                for (int j = left; j < right; j++)
                    neighborBuf.Add(anchorDevs[j]);

                neighborBuf.Sort();
                localBias[i] = neighborBuf[neighborBuf.Count / 2];
            }

            // Compute global sigma of (dev - localBias) across confident targets
            // This measures how tightly peptides cluster around their local center.
            double varSum = 0.0;
            int varCount = 0;
            for (int i = 0; i < targetIndices.Count; i++)
            {
                int idx = targetIndices[i];
                if (!active[idx]) continue;
                var r = results[idx];
                if (r.FdrInfo == null || r.FdrInfo.QValue > 0.01f) continue;
                float dev = r.RtDeviationMinutes;
                float bias = localBias[idx];
                if (float.IsNaN(dev) || float.IsNaN(bias)) continue;
                double residual = dev - bias;
                varSum += residual * residual;
                varCount++;
            }

            sigma = varCount > 1
                ? (float)Math.Sqrt(varSum / varCount)
                : float.NaN;
        }

        private static int LowerBound(float[] arr, float value)
        {
            int lo = 0, hi = arr.Length;
            while (lo < hi)
            {
                int mid = (lo + hi) / 2;
                if (arr[mid] < value) lo = mid + 1; else hi = mid;
            }
            return lo;
        }

        private static int UpperBound(float[] arr, float value)
        {
            int lo = 0, hi = arr.Length;
            while (lo < hi)
            {
                int mid = (lo + hi) / 2;
                if (arr[mid] <= value) lo = mid + 1; else hi = mid;
            }
            return lo;
        }

        /// <summary>
        /// Legacy global stats — kept for reference but no longer used in refinement.
        /// </summary>
        private static void ComputeRtStats(
            List<DiaSearchResult> results,
            List<int> targetIndices,
            bool[] active,
            out float center,
            out float sigma)
        {
            var devs = new List<float>(targetIndices.Count / 4);
            for (int i = 0; i < targetIndices.Count; i++)
            {
                int idx = targetIndices[i];
                if (!active[idx]) continue;
                var r = results[idx];
                if (r.FdrInfo == null || r.FdrInfo.QValue > 0.01f) continue;
                float dev = r.RtDeviationMinutes;
                if (float.IsNaN(dev)) continue;
                devs.Add(dev);
            }

            if (devs.Count < 10) { center = float.NaN; sigma = float.NaN; return; }

            devs.Sort();
            center = devs[devs.Count / 2];

            float mean = 0f;
            for (int i = 0; i < devs.Count; i++) mean += devs[i];
            mean /= devs.Count;

            float variance = 0f;
            for (int i = 0; i < devs.Count; i++)
            {
                float d = devs[i] - mean;
                variance += d * d;
            }
            sigma = MathF.Sqrt(variance / devs.Count);
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
            DiaClassifierType classifierType,
            IDiaClassifier classifier)
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

            // Extract classifier-specific info
            float[] weights = null;
            float bias = 0f;
            float[] importances = null;

            if (classifier is DiaLdaClassifierAdapter ldaAdapter && ldaAdapter.Lda != null)
            {
                weights = (float[])ldaAdapter.Lda.Weights.Clone();
                bias = ldaAdapter.Lda.Bias;
            }

            return new IterationDiagnostics(
                iteration, positiveTrainCount, negativeTrainCount,
                tStats.Mean, tStats.Median, tStats.Q25, tStats.Q75,
                dStats.Mean, dStats.Median, dStats.Q25, dStats.Q75,
                separation, auc, idsAt1Pct,
                weightChange, classifierType, weights, bias, importances);
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
            Console.WriteLine($"  ─── Iteration {d.Iteration} ({d.ClassifierType}) ───");
            Console.WriteLine($"    Training: {d.PositiveTrainingCount:N0} positives, {d.NegativeTrainingCount:N0} negatives");
            Console.WriteLine($"    Target scores:  mean={d.TargetScoreMean:F4}  median={d.TargetScoreMedian:F4}  Q25={d.TargetScoreQ25:F4}  Q75={d.TargetScoreQ75:F4}");
            Console.WriteLine($"    Decoy scores:   mean={d.DecoyScoreMean:F4}  median={d.DecoyScoreMedian:F4}  Q25={d.DecoyScoreQ25:F4}  Q75={d.DecoyScoreQ75:F4}");
            Console.WriteLine($"    Separation:     {d.Separation:F4} (mean_target - mean_decoy)");
            Console.WriteLine($"    AUC:            {d.Auc:F4}");
            Console.WriteLine($"    IDs at 1% FDR:  {d.IdentificationsAt1Pct:N0}");

            if (d.ClassifierType == DiaClassifierType.LinearDiscriminant && d.Weights != null)
            {
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
            else if (d.ClassifierType == DiaClassifierType.GradientBoostedTree)
            {
                Console.WriteLine($"    Classifier:     GBT (200 trees, depth 5, lr 0.05)");
                if (d.FeatureImportances != null)
                {
                    int nFeats = Math.Min(d.FeatureImportances.Length, DiaFeatureVector.ClassifierFeatureCount);
                    var indexed = new (string Name, float Importance)[nFeats];
                    for (int i = 0; i < nFeats; i++)
                        indexed[i] = (DiaFeatureVector.FeatureNames[i], d.FeatureImportances[i]);
                    Array.Sort(indexed, (a, b) => b.Importance.CompareTo(a.Importance));
                    Console.WriteLine("    Feature Importances:");
                    for (int i = 0; i < nFeats; i++)
                        Console.WriteLine($"      {indexed[i].Importance:F4}  {indexed[i].Name}");
                }
            }
        }
    }
}