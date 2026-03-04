// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
// Location: Development/Dia/Phase14GbtSweep.cs

using MassSpectrometry.Dia;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Development.Dia
{
    /// <summary>
    /// Phase 14: GBT Hyperparameter Sweep Utility.
    /// 
    /// Runs LDA baseline + 5 GBT configurations against deep-cloned results.
    /// Each FDR run modifies ClassifierScore and FdrInfo in-place, so results
    /// must be cloned before each config to prevent cross-contamination.
    /// </summary>
    public static class Phase14GbtSweep
    {
        // ════════════════════════════════════════════════════════════════
        //  Sweep Result
        // ════════════════════════════════════════════════════════════════

        /// <summary>Per-config results from the GBT sweep.</summary>
        public readonly struct SweepResult
        {
            public readonly string Name;
            public readonly int IdsAtQ0001;
            public readonly int IdsAtQ0005;
            public readonly int IdsAtQ001;
            public readonly int IdsAtQ005;
            public readonly int IdsAtQ010;
            public readonly int Iterations;
            public readonly long ElapsedMs;

            public SweepResult(string name, int q0001, int q0005, int q001, int q005, int q010,
                int iterations, long elapsedMs)
            {
                Name = name;
                IdsAtQ0001 = q0001;
                IdsAtQ0005 = q0005;
                IdsAtQ001 = q001;
                IdsAtQ005 = q005;
                IdsAtQ010 = q010;
                Iterations = iterations;
                ElapsedMs = elapsedMs;
            }
        }

        // ════════════════════════════════════════════════════════════════
        //  Standard Configurations
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Returns the 5 standard GBT sweep configurations.
        /// Public so the benchmark runner can access config names.
        /// </summary>
        public static DiaFdrEngine.GbtHyperparameters[] GetStandardConfigs()
        {
            return new[]
            {
                new DiaFdrEngine.GbtHyperparameters("A: Conservative",  50, 3, 0.10f, 50, 2.0f, 64, 0.7f),
                new DiaFdrEngine.GbtHyperparameters("B: Moderate",     100, 3, 0.10f, 30, 1.0f, 128, 0.8f),
                new DiaFdrEngine.GbtHyperparameters("C: Deeper",       100, 4, 0.05f, 30, 1.5f, 128, 0.8f),
                new DiaFdrEngine.GbtHyperparameters("D: Shallow-many", 200, 2, 0.10f, 50, 2.0f, 64, 0.7f),
                new DiaFdrEngine.GbtHyperparameters("E: Phase13-orig", 200, 5, 0.05f, 15, 0.5f, 128, 0.85f),
            };
        }

        // ════════════════════════════════════════════════════════════════
        //  Main Sweep Entry Point
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Runs the GBT sweep: 5 configs, each on deep-cloned results.
        /// The original results list is NOT modified.
        /// </summary>
        public static SweepResult[] RunGbtSweep(
            List<DiaSearchResult> originalResults,
            DiaFeatureVector[] features)
        {
            var configs = GetStandardConfigs();
            var sweepResults = new SweepResult[configs.Length];
            float[] thresholds = { 0.001f, 0.005f, 0.01f, 0.05f, 0.10f };

            for (int c = 0; c < configs.Length; c++)
            {
                var config = configs[c];
                Console.WriteLine($"  Config {config.Name}: trees={config.NumTrees} depth={config.MaxDepth} " +
                    $"lr={config.LearningRate:F2} minLeaf={config.MinSamplesLeaf} " +
                    $"L2={config.L2Regularization:F1} subsample={config.SubsampleFraction:F2}");

                // Deep-clone results so FDR doesn't corrupt subsequent runs
                var clonedResults = CloneResults(originalResults);

                var sw = Stopwatch.StartNew();
                var fdrResult = DiaFdrEngine.RunIterativeFdr(clonedResults, features, config);
                sw.Stop();

                // Count IDs at each threshold
                int[] counts = CountIdsAtThresholds(clonedResults, thresholds);

                sweepResults[c] = new SweepResult(
                    config.Name,
                    counts[0], counts[1], counts[2], counts[3], counts[4],
                    fdrResult.IterationsCompleted,
                    sw.ElapsedMilliseconds);

                Console.WriteLine($"    → 1% FDR: {counts[2]:N0} IDs | iterations={fdrResult.IterationsCompleted} | {sw.ElapsedMilliseconds}ms");
            }

            return sweepResults;
        }

        // ════════════════════════════════════════════════════════════════
        //  Summary Table
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Prints side-by-side comparison of LDA baseline vs all GBT configs.
        /// </summary>
        public static void PrintSummaryTable(SweepResult[] sweepResults, int[] ldaCounts)
        {
            Console.WriteLine();
            Console.WriteLine("  GBT Sweep Summary:");
            Console.WriteLine("  ──────────────────────────────────────────────────────────────────────────────");
            Console.WriteLine("  Config              |  q≤0.1%  q≤0.5%  q≤1.0%  q≤5.0%  q≤10%   | Iters | Time");
            Console.WriteLine("  ──────────────────────────────────────────────────────────────────────────────");

            // LDA baseline row
            if (ldaCounts != null && ldaCounts.Length >= 5)
            {
                Console.WriteLine($"  {"LDA baseline",-20} | {ldaCounts[0],7:N0} {ldaCounts[1],7:N0} {ldaCounts[2],7:N0} " +
                    $"{ldaCounts[3],7:N0} {ldaCounts[4],7:N0}  |   —   |   —");
            }

            // GBT rows
            int bestIdx = -1;
            int bestAt1Pct = 0;
            for (int i = 0; i < sweepResults.Length; i++)
            {
                if (sweepResults[i].IdsAtQ001 > bestAt1Pct)
                {
                    bestAt1Pct = sweepResults[i].IdsAtQ001;
                    bestIdx = i;
                }
            }

            for (int i = 0; i < sweepResults.Length; i++)
            {
                var sr = sweepResults[i];
                string marker = (i == bestIdx) ? " ★" : "  ";
                Console.WriteLine($"  {sr.Name,-20}{marker}| {sr.IdsAtQ0001,7:N0} {sr.IdsAtQ0005,7:N0} {sr.IdsAtQ001,7:N0} " +
                    $"{sr.IdsAtQ005,7:N0} {sr.IdsAtQ010,7:N0}  | {sr.Iterations,5} | {sr.ElapsedMs}ms");
            }

            Console.WriteLine("  ──────────────────────────────────────────────────────────────────────────────");

            // Decision recommendation
            if (ldaCounts != null && ldaCounts.Length >= 3)
            {
                int ldaAt1Pct = ldaCounts[2];
                int delta = bestAt1Pct - ldaAt1Pct;

                if (delta > 500)
                    Console.WriteLine($"  ✓ RECOMMENDATION: GBT beats LDA by {delta:N0} IDs. Adopt best config as default.");
                else if (delta > 0)
                    Console.WriteLine($"  ~ GBT marginally better (+{delta:N0}). Consider keeping both; may not justify complexity.");
                else
                    Console.WriteLine($"  ✗ GBT does not beat LDA ({delta:N0}). Non-linear signal may not exist in current features.");
            }
            Console.WriteLine();
        }

        // ════════════════════════════════════════════════════════════════
        //  Deep Clone
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Deep-clones the results list so that in-place FDR modifications
        /// (ClassifierScore, FdrInfo) don't contaminate subsequent sweep runs.
        /// 
        /// Uses the constructor for get-only properties, then copies mutable props.
        /// Feature vectors (DiaFeatureVector[]) are NOT cloned — they are read-only
        /// during FDR and shared across all runs.
        /// </summary>
        private static List<DiaSearchResult> CloneResults(List<DiaSearchResult> originals)
        {
            var cloned = new List<DiaSearchResult>(originals.Count);
            for (int i = 0; i < originals.Count; i++)
            {
                var o = originals[i];
                var c = new DiaSearchResult(
                    o.Sequence, o.ChargeState, o.PrecursorMz, o.WindowId, o.IsDecoy,
                    o.FragmentsQueried, o.LibraryRetentionTime, o.RtWindowStart, o.RtWindowEnd);

                // Copy mutable arrays
                if (o.ExtractedIntensities != null)
                    Array.Copy(o.ExtractedIntensities, c.ExtractedIntensities, o.FragmentsQueried);
                if (o.XicPointCounts != null)
                    Array.Copy(o.XicPointCounts, c.XicPointCounts, o.FragmentsQueried);

                // Copy mutable score/evidence properties
                c.DotProductScore = o.DotProductScore;
                c.SpectralAngleScore = o.SpectralAngleScore;
                c.FragmentsDetected = o.FragmentsDetected;

                // Temporal scoring
                c.ApexScore = o.ApexScore;
                c.TemporalScore = o.TemporalScore;
                c.SpectralAngle = o.SpectralAngle;
                c.MeanFragCorr = o.MeanFragCorr;
                c.MinFragCorr = o.MinFragCorr;
                c.PeakMeanFragCorr = o.PeakMeanFragCorr;
                c.PeakMinFragCorr = o.PeakMinFragCorr;
                c.PeakApexScore = o.PeakApexScore;
                c.PeakWidth = o.PeakWidth;
                c.PeakSymmetry = o.PeakSymmetry;
                c.CandidateCount = o.CandidateCount;
                c.LogTotalIntensity = o.LogTotalIntensity;
                c.IntensityCV = o.IntensityCV;
                c.FragDetRate = o.FragDetRate;
                c.RtDeviationMinutes = o.RtDeviationMinutes;
                c.RtDeviationSquared = o.RtDeviationSquared;

                // Mass accuracy
                c.MeanMassErrorPpm = o.MeanMassErrorPpm;
                c.MedianMassErrorPpm = o.MedianMassErrorPpm;
                c.MassErrorStdPpm = o.MassErrorStdPpm;
                c.MaxAbsMassErrorPpm = o.MaxAbsMassErrorPpm;
                if (o.ApexObservedMzs != null)
                {
                    c.ApexObservedMzs = new float[o.ApexObservedMzs.Length];
                    Array.Copy(o.ApexObservedMzs, c.ApexObservedMzs, o.ApexObservedMzs.Length);
                }

                // Best-fragment reference curve
                c.BestFragCorrelationSum = o.BestFragCorrelationSum;
                c.MedianFragRefCorr = o.MedianFragRefCorr;
                c.MinFragRefCorr = o.MinFragRefCorr;
                c.StdFragRefCorr = o.StdFragRefCorr;

                // Signal ratio
                c.MeanSignalRatioDev = o.MeanSignalRatioDev;
                c.MaxSignalRatioDev = o.MaxSignalRatioDev;
                c.StdSignalRatioDev = o.StdSignalRatioDev;

                // Smoothed correlations & S/N
                c.SmoothedMeanFragCorr = o.SmoothedMeanFragCorr;
                c.SmoothedMinFragCorr = o.SmoothedMinFragCorr;
                c.Log2SignalToNoise = o.Log2SignalToNoise;

                // Migrated features
                c.BestFragWeightedCosine = o.BestFragWeightedCosine;
                c.BestFragIndex = o.BestFragIndex;
                c.BoundarySignalRatio = o.BoundarySignalRatio;
                c.ApexToMeanRatio = o.ApexToMeanRatio;

                // Temporal metadata
                c.TimePointsUsed = o.TimePointsUsed;
                c.ObservedApexRt = o.ObservedApexRt;
                c.DetectedPeakGroup = o.DetectedPeakGroup;

                // Reset FDR fields (these are what we're testing)
                c.ClassifierScore = float.NaN;
                c.FdrInfo = null;

                cloned.Add(c);
            }
            return cloned;
        }

        // ════════════════════════════════════════════════════════════════
        //  Utility
        // ════════════════════════════════════════════════════════════════

        private static int[] CountIdsAtThresholds(List<DiaSearchResult> results, float[] thresholds)
        {
            int[] counts = new int[thresholds.Length];
            for (int i = 0; i < results.Count; i++)
            {
                if (results[i].IsDecoy) continue;
                var fdr = results[i].FdrInfo;
                if (fdr == null) continue;
                for (int t = 0; t < thresholds.Length; t++)
                    if (fdr.QValue <= thresholds[t])
                        counts[t]++;
            }
            return counts;
        }
    }
}