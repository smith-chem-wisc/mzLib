// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// Phase 15, Prompts 2-3 (fixed in Prompt 5): Iterative calibration framework with model selection
// Placement: MassSpectrometry/Dia/Calibration/IterativeRtCalibrator.cs
//
// Prompt 5 fixes:
//   1. Anchor selection uses the same RT coordinate (iRT preferred) as query generation,
//      preventing coordinate mismatch between fitting and prediction.
//   2. Slope-based divergence protection: if slope moves away from 1.0, revert.
//   3. Tighter initial RANSAC inlier threshold (1.0 min instead of 2.0) for bootstrap.
//   4. Console progress logging so the user sees which iteration is running.
//   5. Higher initial ApexScore threshold (0.85) and stricter initial anchor selection
//      to resist contamination from the broad ±5.0 min window.

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using MassSpectrometry.Dia.Calibration;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Implements iterative RT calibration following the progressive refinement strategy
    /// used by DIA-NN and Spectronaut.
    ///
    /// Algorithm:
    /// 1. Ultra-broad initial extraction (±5.0 min or wider)
    /// 2. Select top-K anchors by composite score (robust apex + library cosine)
    /// 3. Fit initial linear model
    /// 4. Re-extract at predicted RT ± Nσ (starting at 3σ)
    /// 5. Re-select anchors with improved confidence
    /// 6. Re-fit model (with automatic model selection: Linear → PiecewiseLinear → LOWESS)
    /// 7. Repeat until convergence (Δσ/σ &lt; threshold) or max iterations
    ///
    /// Model Selection Logic (Prompt 3):
    /// - Iteration 0: always Linear (few anchors, avoid overfitting)
    /// - Iteration 1+: if linear R² &lt; 0.995 and anchors &gt; 200, try PiecewiseLinear
    /// - If PiecewiseLinear shows &gt; 10% σ improvement over Linear, keep it
    /// - LOWESS tried only if PiecewiseLinear R² &lt; 0.990 (extreme non-linearity)
    ///
    /// Prompt 5 critical fixes:
    /// - Anchor selection resolves the same RT coordinate as query generation
    ///   (iRT preferred over RetentionTime) to prevent coordinate mismatch
    /// - Slope divergence protection: reverts if slope moves away from 1.0
    /// - Tighter bootstrap RANSAC (1.0 min) to reject contaminated anchors
    /// </summary>
    public class IterativeRtCalibrator
    {
        // ── Configuration ──────────────────────────────────────────────────

        /// <summary>Maximum iterations including bootstrap. Default 4.</summary>
        public int MaxIterations { get; set; } = 4;

        /// <summary>Stop when Δσ/σ &lt; this threshold (5%). Default 0.05.</summary>
        public double ConvergenceThreshold { get; set; } = 0.05;

        /// <summary>Window = predicted ± N × σ. Default 3.0.</summary>
        public double SigmaMultiplier { get; set; } = 3.0;

        /// <summary>Maximum anchors for iteration 0 (bootstrap). Default 500.</summary>
        public int InitialTopK { get; set; } = 500;

        /// <summary>Strict ApexScore threshold for bootstrap iteration. Default 0.85.</summary>
        public float InitialApexScoreThreshold { get; set; } = 0.85f;

        /// <summary>Relaxed ApexScore threshold for later iterations. Default 0.5.</summary>
        public float RefinedApexScoreThreshold { get; set; } = 0.5f;

        /// <summary>Floor on window half-width to prevent excessively tight windows. Default 0.3 min.</summary>
        public double MinWindowHalfWidthMinutes { get; set; } = 0.3;

        /// <summary>Minimum anchors required to proceed with fitting. Default 20.</summary>
        public int MinAnchorCount { get; set; } = 20;

        /// <summary>
        /// Whether to enable automatic non-linear model selection.
        /// When false, only linear models are used (Prompt 2 behavior).
        /// When true, PiecewiseLinear and LOWESS models are tried when appropriate.
        /// Default true (Prompt 3).
        /// </summary>
        public bool EnableNonLinearModelSelection { get; set; } = true;

        /// <summary>
        /// R² threshold below which PiecewiseLinear is attempted on iterations 1+.
        /// Default 0.995.
        /// </summary>
        public double PiecewiseLinearRSquaredThreshold { get; set; } = 0.995;

        /// <summary>
        /// R² threshold below which LOWESS is attempted. Default 0.990.
        /// </summary>
        public double LowessRSquaredThreshold { get; set; } = 0.990;

        /// <summary>
        /// Minimum σ improvement (fraction) for a non-linear model to be preferred.
        /// Default 0.10 (10%).
        /// </summary>
        public double NonLinearSigmaImprovementThreshold { get; set; } = 0.10;

        /// <summary>Number of segments for PiecewiseLinear model. Default 8.</summary>
        public int PiecewiseLinearSegments { get; set; } = 8;

        /// <summary>Bandwidth parameter for LOWESS fitting. Default 0.3.</summary>
        public double LowessBandwidth { get; set; } = 0.3;

        // ── Main Entry Point ───────────────────────────────────────────────

        /// <summary>
        /// Runs the full iterative calibration loop.
        /// Returns the final calibration model (as the sealed production type),
        /// the internal model wrapper (for local σ access),
        /// the per-iteration log, and the last-iteration extraction results.
        /// </summary>
        public (RtCalibrationModel Model, IRtCalibrationModel DetailedModel,
                List<CalibrationIterationLog> Log, List<DiaSearchResult> Results) Calibrate(
            IList<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            DiaSearchParameters baseParameters,
            DiaExtractionOrchestrator orchestrator)
        {
            var log = new List<CalibrationIterationLog>();
            IRtCalibrationModel currentModel = null;
            List<DiaSearchResult> currentResults = null;
            double previousSigma = double.MaxValue;
            double previousSlope = double.NaN;

            // Build a lookup from result index to precursor index for RT coordinate resolution.
            // This is built once and reused across iterations.
            // (Not needed here — we build it inside SelectAnchors using PrecursorQueryGroup mapping.)

            // ── Bootstrap (Iteration 0): Progressive Widening ──────────────
            // Instead of starting with the full broad window (e.g. ±5.0 min) which
            // produces heavily contaminated anchors, start narrow (±1.0 min) and
            // widen only if we can't find enough high-quality anchors. At ±1.0 min,
            // true matches have their apex within the window, giving clean anchors
            // for the initial model fit. This eliminates the systematic contamination
            // that caused slope=0.85 with the broad-window bootstrap.
            double[] bootstrapWindows = { 1.0, 1.5, 2.0, 3.0, 5.0 };

            for (int iteration = 0; iteration < MaxIterations; iteration++)
            {
                Console.WriteLine($"  [Calibration] Iteration {iteration} starting...");
                var extractSw = new Stopwatch();
                var fitSw = new Stopwatch();
                var selectSw = new Stopwatch();

                DiaLibraryQueryGenerator.GenerationResult genResult = default;
                ExtractionResult extractionResult = null;
                List<DiaSearchResult> results = null;
                List<(double LibraryRt, double ObservedRt, double QualityScore)> anchors = null;

                if (iteration == 0)
                {
                    // ── Progressive widening bootstrap ────────────────────────
                    // Try increasingly wide windows until we get enough anchors.
                    foreach (double tryWindow in bootstrapWindows)
                    {
                        extractSw.Reset();
                        selectSw.Reset();

                        Console.WriteLine($"  [Calibration] Bootstrap: trying ±{tryWindow:F1} min window...");
                        extractSw.Start();

                        // Create modified parameters with this window width
                        var bootstrapParams = new DiaSearchParameters
                        {
                            PpmTolerance = baseParameters.PpmTolerance,
                            RtToleranceMinutes = (float)tryWindow,
                            MinFragmentsRequired = baseParameters.MinFragmentsRequired,
                            MinScoreThreshold = baseParameters.MinScoreThreshold,
                            MaxThreads = baseParameters.MaxThreads,
                            PreferGpu = baseParameters.PreferGpu,
                            CalibratedWindowSigmaMultiplier = baseParameters.CalibratedWindowSigmaMultiplier,
                            ScoringStrategy = baseParameters.ScoringStrategy,
                            NonlinearPower = baseParameters.NonlinearPower
                        };

                        genResult = DiaLibraryQueryGenerator.Generate(precursors, scanIndex, bootstrapParams);
                        extractionResult = orchestrator.ExtractAll(genResult.Queries);

                        results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                            precursors, genResult, extractionResult, baseParameters, scanIndex);
                        extractSw.Stop();

                        Console.WriteLine($"  [Calibration] Bootstrap ±{tryWindow:F1} min: extraction={extractSw.ElapsedMilliseconds}ms, results={results.Count:N0}");

                        selectSw.Start();
                        anchors = SelectAnchors(
                            results, precursors, genResult.PrecursorGroups,
                            InitialApexScoreThreshold, InitialTopK,
                            requireMinFragDetRate: false);
                        selectSw.Stop();

                        Console.WriteLine($"  [Calibration] Bootstrap ±{tryWindow:F1} min: {anchors.Count} anchors (threshold={InitialApexScoreThreshold:F2})");

                        if (anchors.Count >= InitialTopK || anchors.Count >= MinAnchorCount * 5)
                            break; // Got enough good anchors, no need to widen
                    }
                }
                else
                {
                    // ── Refinement: calibrated extraction ──────────────────────
                    extractSw.Start();

                    double windowHalfWidth = Math.Max(
                        currentModel.SigmaMinutes * SigmaMultiplier,
                        MinWindowHalfWidthMinutes);

                    genResult = GenerateCalibratedWithExplicitWindow(
                        precursors, scanIndex, baseParameters, currentModel, windowHalfWidth);
                    extractionResult = orchestrator.ExtractAll(genResult.Queries);

                    results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                        precursors, genResult, extractionResult, baseParameters, scanIndex);
                    extractSw.Stop();

                    Console.WriteLine($"  [Calibration] Iteration {iteration}: extraction={extractSw.ElapsedMilliseconds}ms, results={results.Count:N0}");

                    selectSw.Start();
                    anchors = SelectAnchors(
                        results, precursors, genResult.PrecursorGroups,
                        RefinedApexScoreThreshold, int.MaxValue,
                        requireMinFragDetRate: true);
                    selectSw.Stop();

                    Console.WriteLine($"  [Calibration] Iteration {iteration}: {anchors.Count} anchors (threshold={RefinedApexScoreThreshold:F2})");
                }

                // Safety: not enough anchors
                if (anchors.Count < MinAnchorCount)
                {
                    log.Add(CreateLogEntry(iteration, anchors.Count, currentModel,
                        currentModel != null
                            ? Math.Max(currentModel.SigmaMinutes * SigmaMultiplier, MinWindowHalfWidthMinutes)
                            : baseParameters.RtToleranceMinutes,
                        extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed,
                        "InsufficientAnchors"));

                    if (currentModel != null)
                        break;
                    else
                        return (null, null, log, results);
                }

                // ── Model Fitting ──────────────────────────────────────────
                fitSw.Start();
                var libraryRts = anchors.Select(a => a.LibraryRt).ToArray();
                var observedRts = anchors.Select(a => a.ObservedRt).ToArray();

                // RANSAC inlier threshold — looser for bootstrap since narrow window
                // already provides clean anchors
                double inlierThreshold = (iteration == 0)
                    ? 1.0
                    : Math.Max(1.5 * previousSigma, 0.3);

                // Fit linear model (always — baseline and for iteration 0)
                var linearModel = FitLinearWithRansac(libraryRts, observedRts, inlierThreshold);

                if (linearModel == null || !linearModel.IsReliable)
                {
                    // Try with looser threshold before giving up
                    if (iteration == 0)
                    {
                        Console.WriteLine($"  [Calibration] Iteration {iteration}: tight RANSAC failed, retrying with 2.0 min threshold...");
                        linearModel = FitLinearWithRansac(libraryRts, observedRts, 2.0);
                    }

                    if (linearModel == null || !linearModel.IsReliable)
                    {
                        log.Add(CreateLogEntry(iteration, anchors.Count, currentModel,
                            currentModel != null
                                ? Math.Max(currentModel.SigmaMinutes * SigmaMultiplier, MinWindowHalfWidthMinutes)
                                : baseParameters.RtToleranceMinutes,
                            extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed,
                            "FitFailed_Reverted"));

                        if (currentModel != null)
                            break;
                        else
                            return (null, null, log, results);
                    }
                }

                // Wrap the linear model
                IRtCalibrationModel bestModel = new LinearRtModelWrapper(linearModel);
                string modelTypeLabel = "Linear";

                // ── Model Selection (Prompt 3) ─────────────────────────────
                if (EnableNonLinearModelSelection && iteration > 0 && anchors.Count > 200)
                {
                    bestModel = SelectBestModel(
                        libraryRts, observedRts, linearModel,
                        out modelTypeLabel);
                }

                fitSw.Stop();

                // ── Convergence / Divergence Check ─────────────────────────
                double newSigma = bestModel.SigmaMinutes;
                double newSlope = bestModel.Slope;

                Console.WriteLine($"  [Calibration] Iteration {iteration}: slope={newSlope:F4}, σ={newSigma:F3}, R²={bestModel.RSquared:F4}, model={modelTypeLabel}");

                // Divergence protection: if σ increases by >10%, revert
                if (iteration > 0 && newSigma > previousSigma * 1.1)
                {
                    Console.WriteLine($"  [Calibration] Iteration {iteration}: σ diverged ({newSigma:F3} > {previousSigma * 1.1:F3}), reverting to previous model.");
                    log.Add(CreateLogEntry(iteration, anchors.Count, bestModel,
                        Math.Max(newSigma * SigmaMultiplier, MinWindowHalfWidthMinutes),
                        extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed,
                        modelTypeLabel + "_Diverged_Reverted"));
                    break;
                }

                // PROMPT 5 FIX: Slope divergence protection.
                // For same-instrument libraries, slope should be near 1.0.
                // If slope is moving AWAY from 1.0 across iterations, something is wrong.
                if (iteration > 0 && !double.IsNaN(previousSlope))
                {
                    double currentDistFromOne = Math.Abs(newSlope - 1.0);
                    double previousDistFromOne = Math.Abs(previousSlope - 1.0);

                    if (currentDistFromOne > previousDistFromOne + 0.02)
                    {
                        Console.WriteLine($"  [Calibration] Iteration {iteration}: slope diverging from 1.0 " +
                                          $"(|{newSlope:F4}-1.0|={currentDistFromOne:F3} > |{previousSlope:F4}-1.0|={previousDistFromOne:F3}), " +
                                          $"reverting to previous model.");
                        log.Add(CreateLogEntry(iteration, anchors.Count, bestModel,
                            Math.Max(newSigma * SigmaMultiplier, MinWindowHalfWidthMinutes),
                            extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed,
                            modelTypeLabel + "_SlopeDiverged_Reverted"));
                        break;
                    }
                }

                // Accept the new model
                currentModel = bestModel;
                currentResults = results;

                double windowHW = Math.Max(newSigma * SigmaMultiplier, MinWindowHalfWidthMinutes);
                log.Add(CreateLogEntry(iteration, anchors.Count, bestModel, windowHW,
                    extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, modelTypeLabel));

                // Convergence check: |σ_new - σ_old| / σ_old < threshold
                if (iteration > 0 && previousSigma > 0 && !double.IsInfinity(previousSigma))
                {
                    double relativeDelta = Math.Abs(newSigma - previousSigma) / previousSigma;
                    if (relativeDelta < ConvergenceThreshold)
                    {
                        Console.WriteLine($"  [Calibration] Converged at iteration {iteration} (Δσ/σ = {relativeDelta:F4} < {ConvergenceThreshold})");
                        break;
                    }
                }

                previousSigma = newSigma;
                previousSlope = newSlope;
            }

            if (currentModel == null)
                return (null, null, log, currentResults);

            return (currentModel.ToRtCalibrationModel(), currentModel, log, currentResults);
        }

        // ── Model Selection (Prompt 3, Task 4) ────────────────────────────

        /// <summary>
        /// Selects the best calibration model type based on the linear model's R² and
        /// the improvement offered by non-linear alternatives.
        /// </summary>
        private IRtCalibrationModel SelectBestModel(
            double[] libraryRts,
            double[] observedRts,
            RtCalibrationModel linearModel,
            out string selectedModelType)
        {
            selectedModelType = "Linear";
            IRtCalibrationModel bestModel = new LinearRtModelWrapper(linearModel);
            double bestSigma = linearModel.SigmaMinutes;

            // If linear R² is already very good, no need for non-linear
            if (linearModel.RSquared >= PiecewiseLinearRSquaredThreshold)
                return bestModel;

            try
            {
                // ── Try PiecewiseLinear ────────────────────────────────────
                var pwlModel = PiecewiseLinearRtModel.Fit(libraryRts, observedRts, PiecewiseLinearSegments);

                if (pwlModel != null && pwlModel.IsReliable)
                {
                    double sigmaImprovement = (bestSigma - pwlModel.SigmaMinutes) / bestSigma;
                    if (sigmaImprovement > NonLinearSigmaImprovementThreshold)
                    {
                        bestModel = pwlModel;
                        bestSigma = pwlModel.SigmaMinutes;
                        selectedModelType = "PiecewiseLinear";
                    }
                }

                // ── Try LOWESS if still insufficient ──────────────────────
                if (bestModel.RSquared < LowessRSquaredThreshold)
                {
                    var lowessModel = LowessRtModel.Fit(
                        libraryRts, observedRts, LowessBandwidth, enforceMonotonic: true);

                    if (lowessModel != null && lowessModel.IsReliable)
                    {
                        double sigmaImpLowess = (bestSigma - lowessModel.SigmaMinutes) / bestSigma;
                        if (sigmaImpLowess > NonLinearSigmaImprovementThreshold)
                        {
                            bestModel = lowessModel;
                            bestSigma = lowessModel.SigmaMinutes;
                            selectedModelType = "Lowess";
                        }
                    }
                }
            }
            catch
            {
                // Non-linear fitting failures should never crash calibration
            }

            return bestModel;
        }

        // ── Anchor Selection ───────────────────────────────────────────────

        /// <summary>
        /// Selects high-quality calibration anchors from extraction results.
        ///
        /// CRITICAL (Prompt 5 fix): This overload uses the precursors list to resolve
        /// the library RT coordinate that the calibration model will be applied to.
        /// It uses the SAME coordinate as query generation (iRT preferred, then
        /// RetentionTime) to prevent coordinate mismatch between fitting and prediction.
        ///
        /// The mapping from result → precursor uses Sequence+ChargeState matching,
        /// since results are not guaranteed to be parallel to PrecursorQueryGroups
        /// (some groups may be skipped during assembly).
        ///
        /// Eligibility: target peptides only, valid library RT and observed apex RT,
        /// ApexScore ≥ minApexScore, and (optionally) FragDetRate ≥ 0.5.
        ///
        /// Composite quality:
        ///   QualityScore = ApexScore × FragDetRate × MassAccuracyTerm
        ///   where MassAccuracyTerm = 1 / (1 + |MeanMassErrorPpm| / 10)
        ///   (NaN mass error → neutral term of 1.0)
        ///
        /// RT-balanced selection when maxAnchors limits apply:
        /// - Divide library RT range into 10 equal bins
        /// - Reserve at least maxAnchors/20 slots per bin
        /// - Fill remaining with highest-quality unused anchors
        /// </summary>
        public static List<(double LibraryRt, double ObservedRt, double QualityScore)> SelectAnchors(
            List<DiaSearchResult> results,
            IList<LibraryPrecursorInput> precursors,
            DiaLibraryQueryGenerator.PrecursorQueryGroup[] groups,
            float minApexScore,
            int maxAnchors,
            bool requireMinFragDetRate = true)
        {
            // Build a lookup from (Sequence, Charge) → precursor for RT coordinate resolution.
            // This handles the fact that results may not be parallel to groups.
            Dictionary<(string Seq, int Charge), LibraryPrecursorInput> precursorLookup = null;
            if (precursors != null && precursors.Count > 0)
            {
                precursorLookup = new Dictionary<(string, int), LibraryPrecursorInput>(precursors.Count);
                for (int i = 0; i < precursors.Count; i++)
                {
                    var p = precursors[i];
                    if (p.IsDecoy) continue; // Only need targets for anchors
                    var key = (p.Sequence, p.ChargeState);
                    if (!precursorLookup.ContainsKey(key))
                        precursorLookup[key] = p;
                }
            }

            var candidates = new List<(double LibraryRt, double ObservedRt, double QualityScore)>();

            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];

                // Target peptides only
                if (r.IsDecoy) continue;

                // Valid observed apex RT
                if (float.IsNaN(r.ObservedApexRt)) continue;
                double obsRt = r.ObservedApexRt;

                // ApexScore threshold
                if (float.IsNaN(r.ApexScore) || r.ApexScore < minApexScore) continue;

                // FragDetRate threshold
                float fragDetRate = r.FragmentDetectionRate;
                if (requireMinFragDetRate && fragDetRate < 0.5f) continue;

                // ── Resolve library RT coordinate ──────────────────────
                // CRITICAL: must match what GenerateCalibratedWithExplicitWindow uses.
                // That method uses p.IrtValue (preferred) then p.RetentionTime.
                double libRt;
                if (precursorLookup != null &&
                    precursorLookup.TryGetValue((r.Sequence, r.ChargeState), out var precursor))
                {
                    if (precursor.IrtValue.HasValue && !double.IsNaN(precursor.IrtValue.Value))
                        libRt = precursor.IrtValue.Value;
                    else if (precursor.RetentionTime.HasValue && !double.IsNaN(precursor.RetentionTime.Value))
                        libRt = precursor.RetentionTime.Value;
                    else
                        continue; // No usable RT coordinate
                }
                else
                {
                    // Fallback: use result's LibraryRetentionTime (old behavior)
                    if (!r.LibraryRetentionTime.HasValue || double.IsNaN(r.LibraryRetentionTime.Value))
                        continue;
                    libRt = r.LibraryRetentionTime.Value;
                }

                // Composite quality score
                double massAccuracyTerm = 1.0;
                if (!float.IsNaN(r.MeanMassErrorPpm))
                    massAccuracyTerm = 1.0 / (1.0 + Math.Abs(r.MeanMassErrorPpm) / 10.0);

                double quality = r.ApexScore * fragDetRate * massAccuracyTerm;
                candidates.Add((libRt, obsRt, quality));
            }

            if (candidates.Count == 0)
                return candidates;

            if (candidates.Count <= maxAnchors)
            {
                candidates.Sort((a, b) => b.QualityScore.CompareTo(a.QualityScore));
                return candidates;
            }

            return SelectWithRtBalance(candidates, maxAnchors, numBins: 10);
        }

        /// <summary>
        /// Legacy overload without groups — uses r.LibraryRetentionTime directly.
        /// Retained for backward compatibility with callers that don't have groups.
        /// </summary>
        public static List<(double LibraryRt, double ObservedRt, double QualityScore)> SelectAnchors(
            List<DiaSearchResult> results,
            IList<LibraryPrecursorInput> precursors,
            float minApexScore,
            int maxAnchors,
            bool requireMinFragDetRate = true)
        {
            return SelectAnchors(results, precursors, null, minApexScore, maxAnchors, requireMinFragDetRate);
        }

        private static List<(double LibraryRt, double ObservedRt, double QualityScore)> SelectWithRtBalance(
            List<(double LibraryRt, double ObservedRt, double QualityScore)> candidates,
            int maxAnchors,
            int numBins)
        {
            double minRt = double.MaxValue, maxRt = double.MinValue;
            foreach (var c in candidates)
            {
                if (c.LibraryRt < minRt) minRt = c.LibraryRt;
                if (c.LibraryRt > maxRt) maxRt = c.LibraryRt;
            }

            double binWidth = (maxRt - minRt) / numBins;
            if (binWidth < 1e-6) binWidth = 1e-6;

            var bins = new List<(double LibraryRt, double ObservedRt, double QualityScore)>[numBins];
            for (int b = 0; b < numBins; b++)
                bins[b] = new List<(double, double, double)>();

            foreach (var c in candidates)
            {
                int b = Math.Min((int)((c.LibraryRt - minRt) / binWidth), numBins - 1);
                bins[b].Add(c);
            }

            for (int b = 0; b < numBins; b++)
                bins[b].Sort((a, x) => x.QualityScore.CompareTo(a.QualityScore));

            int perBinReserve = maxAnchors / (numBins * 2);
            var selected = new List<(double LibraryRt, double ObservedRt, double QualityScore)>(maxAnchors);
            var usedCounts = new int[numBins];

            for (int b = 0; b < numBins; b++)
            {
                int take = Math.Min(perBinReserve, bins[b].Count);
                for (int i = 0; i < take; i++)
                    selected.Add(bins[b][i]);
                usedCounts[b] = take;
            }

            int remaining = maxAnchors - selected.Count;
            if (remaining > 0)
            {
                var unused = new List<(double LibraryRt, double ObservedRt, double QualityScore)>();
                for (int b = 0; b < numBins; b++)
                    for (int i = usedCounts[b]; i < bins[b].Count; i++)
                        unused.Add(bins[b][i]);

                unused.Sort((a, x) => x.QualityScore.CompareTo(a.QualityScore));
                int take = Math.Min(remaining, unused.Count);
                for (int i = 0; i < take; i++)
                    selected.Add(unused[i]);
            }

            return selected;
        }

        // ── Linear Model Fitting with RANSAC ───────────────────────────────

        private static RtCalibrationModel FitLinearWithRansac(
            double[] libraryRts, double[] observedRts, double inlierThreshold)
        {
            var fitOptions = new RtCalibrationFitter.FitOptions
            {
                UseRansac = true,
                RansacInlierThresholdMinutes = inlierThreshold,
                RansacMinInlierFraction = 0.5
            };

            // Use IList<double> overload to avoid ambiguity with ReadOnlySpan
            return RtCalibrationFitter.Fit(
                (IList<double>)libraryRts,
                (IList<double>)observedRts,
                fitOptions);
        }

        // ── Calibrated Query Generation ────────────────────────────────────

        /// <summary>
        /// Generates calibrated queries with an explicit window half-width, bypassing
        /// the CalibratedWindowSigmaMultiplier in DiaSearchParameters. This is needed
        /// during iterative calibration where the window is computed from the current
        /// model's σ and the calibrator's SigmaMultiplier.
        ///
        /// Follows the same two-pass (count, then fill) pattern as
        /// DiaLibraryQueryGenerator.GenerateCalibrated().
        ///
        /// RT coordinate: uses iRT (preferred) then RetentionTime, matching
        /// the anchor selection logic.
        /// </summary>
        private static DiaLibraryQueryGenerator.GenerationResult GenerateCalibratedWithExplicitWindow(
            IList<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            DiaSearchParameters parameters,
            IRtCalibrationModel model,
            double windowHalfWidthMinutes)
        {
            float ppmTolerance = parameters.PpmTolerance;
            float globalRtMin = scanIndex.GetGlobalRtMin();
            float globalRtMax = scanIndex.GetGlobalRtMax();

            // First pass: count
            int totalQueryCount = 0;
            int skippedNoWindow = 0;
            int skippedNoFragments = 0;

            for (int i = 0; i < precursors.Count; i++)
            {
                var p = precursors[i];
                int windowId = scanIndex.FindWindowForPrecursorMz(p.PrecursorMz);
                if (windowId < 0) { skippedNoWindow++; continue; }
                if (p.FragmentCount == 0) { skippedNoFragments++; continue; }
                totalQueryCount += p.FragmentCount;
            }

            var queries = new FragmentQuery[totalQueryCount];
            var groups = new List<DiaLibraryQueryGenerator.PrecursorQueryGroup>(
                precursors.Count - skippedNoWindow - skippedNoFragments);

            // Second pass: fill
            int queryIndex = 0;
            for (int i = 0; i < precursors.Count; i++)
            {
                var p = precursors[i];
                int windowId = scanIndex.FindWindowForPrecursorMz(p.PrecursorMz);
                if (windowId < 0 || p.FragmentCount == 0)
                    continue;

                float rtMin, rtMax;

                // Use iRT if available, else library RT, else full range
                if (p.IrtValue.HasValue)
                {
                    float predictedRt = (float)model.ToMinutes(p.IrtValue.Value);
                    rtMin = (float)(predictedRt - windowHalfWidthMinutes);
                    rtMax = (float)(predictedRt + windowHalfWidthMinutes);
                }
                else if (p.RetentionTime.HasValue)
                {
                    float predictedRt = (float)model.ToMinutes(p.RetentionTime.Value);
                    rtMin = (float)(predictedRt - windowHalfWidthMinutes);
                    rtMax = (float)(predictedRt + windowHalfWidthMinutes);
                }
                else
                {
                    rtMin = globalRtMin;
                    rtMax = globalRtMax;
                }

                int groupOffset = queryIndex;

                for (int f = 0; f < p.FragmentCount; f++)
                {
                    queries[queryIndex] = new FragmentQuery(
                        targetMz: p.FragmentMzs[f],
                        tolerancePpm: ppmTolerance,
                        rtMin: rtMin,
                        rtMax: rtMax,
                        windowId: windowId,
                        queryId: queryIndex);
                    queryIndex++;
                }

                groups.Add(new DiaLibraryQueryGenerator.PrecursorQueryGroup(
                    inputIndex: i,
                    queryOffset: groupOffset,
                    queryCount: p.FragmentCount,
                    windowId: windowId,
                    rtMin: rtMin,
                    rtMax: rtMax));
            }

            return new DiaLibraryQueryGenerator.GenerationResult(
                queries, groups.ToArray(), skippedNoWindow, skippedNoFragments);
        }

        // ── Logging Helpers ────────────────────────────────────────────────

        private static CalibrationIterationLog CreateLogEntry(
            int iteration,
            int anchorCount,
            IRtCalibrationModel model,
            double windowHalfWidth,
            TimeSpan extractionTime,
            TimeSpan fitTime,
            TimeSpan selectionTime,
            string modelType)
        {
            return new CalibrationIterationLog
            {
                Iteration = iteration,
                AnchorCount = anchorCount,
                Slope = model?.Slope ?? double.NaN,
                Intercept = model?.Intercept ?? double.NaN,
                SigmaMinutes = model?.SigmaMinutes ?? double.NaN,
                RSquared = model?.RSquared ?? double.NaN,
                WindowHalfWidthMinutes = windowHalfWidth,
                ModelType = modelType,
                ExtractionTime = extractionTime,
                FitTime = fitTime,
                SelectionTime = selectionTime
            };
        }
    }

    // ── Per-Iteration Diagnostic Log ───────────────────────────────────────

    /// <summary>
    /// Captures diagnostics for a single iteration of the calibration loop.
    /// Used by the benchmark runner (Prompt 5) to produce the convergence table.
    /// </summary>
    public class CalibrationIterationLog
    {
        public int Iteration { get; set; }
        public int AnchorCount { get; set; }
        public double Slope { get; set; }
        public double Intercept { get; set; }
        public double SigmaMinutes { get; set; }
        public double RSquared { get; set; }
        public double WindowHalfWidthMinutes { get; set; }
        public string ModelType { get; set; }
        public TimeSpan ExtractionTime { get; set; }
        public TimeSpan FitTime { get; set; }
        public TimeSpan SelectionTime { get; set; }
        public TimeSpan TotalTime => ExtractionTime + FitTime + SelectionTime;

        public override string ToString()
        {
            return $"Iter {Iteration}: {AnchorCount} anchors, slope={Slope:F4}, σ={SigmaMinutes:F3} min, " +
                   $"R²={RSquared:F4}, window=±{WindowHalfWidthMinutes:F2} min, model={ModelType}";
        }
    }
}