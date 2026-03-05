// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// Phase 15, Prompts 2-3 (fixed in Prompt 5): Iterative calibration framework with model selection
// Phase 21: Offset-detection bootstrap — finds the bulk iRT→RT shift before narrowing.
// Placement: MassSpectrometry/Dia/Calibration/IterativeRtCalibrator.cs

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
    /// 1. Offset-detection bootstrap:
    ///    a. Extract a random subsample of precursors using the full run RT range
    ///    b. Score results, compute per-result offset = ObservedApexRt - LibraryRt
    ///    c. Find the mode of the offset distribution via kernel density
    ///    d. Re-extract ALL precursors at LibraryRt + modeOffset ± bootstrapWindow
    ///    e. Fit initial linear model from anchors
    /// 2. Re-extract at predicted RT ± Nσ (starting at 3σ)
    /// 3. Re-select anchors with improved confidence
    /// 4. Re-fit model (with automatic model selection: Linear → PiecewiseLinear → LOWESS)
    /// 5. Repeat until convergence (Δσ/σ &lt; threshold) or max iterations
    ///
    /// The offset-detection step solves the fundamental problem where iRT values are in
    /// arbitrary units (e.g. -30 to +150) while experimental RT is 0–90 min. The old
    /// progressive-widening bootstrap (±1, ±1.5, ±2, ±3, ±5 min) assumes offset ≈ 0,
    /// which fails completely when the true offset is tens of minutes.
    ///
    /// Model Selection Logic:
    /// - Iteration 0: always Linear (few anchors, avoid overfitting)
    /// - Iteration 1+: if linear R² &lt; 0.995 and anchors &gt; 200, try PiecewiseLinear
    /// - If PiecewiseLinear shows &gt; 10% σ improvement over Linear, keep it
    /// - LOWESS tried only if PiecewiseLinear R² &lt; 0.990 (extreme non-linearity)
    /// </summary>
    public class IterativeRtCalibrator
    {
        // ── Configuration ──────────────────────────────────────────────────

        /// <summary>Maximum iterations including bootstrap. Default 4.</summary>
        public int MaxIterations { get; set; } = 6;

        /// <summary>Stop when Δσ/σ &lt; this threshold (5%). Default 0.05.</summary>
        public double ConvergenceThreshold { get; set; } = 0.02;

        /// <summary>Window = predicted ± N × σ. Default 3.0.</summary>
        public double SigmaMultiplier { get; set; } = 4.0;

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
        /// Number of precursors to subsample for offset detection.
        /// A random subset is used for the global scan to keep extraction fast.
        /// Default 5000. Set to int.MaxValue to use all precursors.
        /// </summary>
        public int OffsetDetectionSubsampleSize { get; set; } = 5000;

        /// <summary>
        /// Half-width of the extraction window used after offset detection.
        /// This window is applied as LibraryRt + modeOffset ± BootstrapWindowMinutes.
        /// Should be wide enough to catch the true apex despite residual spread,
        /// but narrow enough to exclude most false matches.
        /// Default 2.0 min.
        /// </summary>
        public double BootstrapWindowMinutes { get; set; } = 2.0;

        /// <summary>
        /// Minimum number of high-quality results required from the offset-detection
        /// subsample to trust the detected mode. If fewer results are found, falls
        /// back to the legacy progressive-widening bootstrap.
        /// Default 50.
        /// </summary>
        public int OffsetDetectionMinResults { get; set; } = 50;

        /// <summary>
        /// Bandwidth for kernel density estimation of the offset distribution, in minutes.
        /// Controls how sharply the mode is detected. Default 0.5 min.
        /// </summary>
        public double OffsetKdeBandwidthMinutes { get; set; } = 0.5;

        /// <summary>
        /// Whether to enable automatic non-linear model selection.
        /// When false, only linear models are used.
        /// Default true.
        /// </summary>
        public bool EnableNonLinearModelSelection { get; set; } = true;

        /// <summary>R² threshold below which PiecewiseLinear is attempted. Default 0.995.</summary>
        public double PiecewiseLinearRSquaredThreshold { get; set; } = 0.995;

        /// <summary>R² threshold below which LOWESS is attempted. Default 0.990.</summary>
        public double LowessRSquaredThreshold { get; set; } = 0.990;

        /// <summary>Minimum σ improvement (fraction) for a non-linear model to be preferred. Default 0.10.</summary>
        public double NonLinearSigmaImprovementThreshold { get; set; } = 0.10;

        /// <summary>Number of segments for PiecewiseLinear model. Default 8.</summary>
        public int PiecewiseLinearSegments { get; set; } = 8;

        /// <summary>Bandwidth parameter for LOWESS fitting. Default 0.3.</summary>
        public double LowessBandwidth { get; set; } = 0.3;

        // ── Main Entry Point ───────────────────────────────────────────────

        /// <summary>
        /// Runs the full iterative calibration loop.
        /// Returns the final calibration model, the internal model wrapper,
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

            // ── Bootstrap (Iteration 0): Offset Detection ─────────────────
            // Phase 21: detect the bulk iRT→RT offset before attempting any
            // narrow-window extraction.
            DiaLibraryQueryGenerator.GenerationResult genResult = default;
            ExtractionResult extractionResult = null;
            List<DiaSearchResult> bootstrapResults = null;

            {
                var extractSw = new Stopwatch();
                var fitSw = new Stopwatch();
                var selectSw = new Stopwatch();

                Console.WriteLine("  [Calibration] Iteration 0 starting...");
                Console.WriteLine("  [Calibration] Bootstrap: offset-detection scan...");

                extractSw.Start();

                // Step A: Subsample precursors for global scan
                var subsample = BuildSubsample(precursors, OffsetDetectionSubsampleSize);

                // Step B: Extract subsample over full RT range
                var subsampleParams = CloneWithRtTolerance(baseParameters,
                    (scanIndex.GetGlobalRtMax() - scanIndex.GetGlobalRtMin()) / 2f + 1f);

                var subGen = DiaLibraryQueryGenerator.Generate(subsample, scanIndex, subsampleParams);
                var subExtract = orchestrator.ExtractAll(subGen.Queries);
                var subResults = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    subsample, subGen, subExtract, baseParameters, scanIndex);

                extractSw.Stop();
                Console.WriteLine($"  [Calibration] Offset scan: {subResults.Count} results from {subsample.Count} precursors ({extractSw.ElapsedMilliseconds}ms)");

                // Step C: Detect mode of offset distribution
                double modeOffset = double.NaN;
                bool offsetDetected = false;

                var highQualitySub = subResults
                    .Where(r => !r.IsDecoy && !float.IsNaN(r.ApexScore) && r.ApexScore >= InitialApexScoreThreshold
                                && !float.IsNaN(r.ObservedApexRt))
                    .ToList();

                Console.WriteLine($"  [Calibration] Offset scan: {highQualitySub.Count} high-quality results (ApexScore≥{InitialApexScoreThreshold:F2})");

                if (highQualitySub.Count >= OffsetDetectionMinResults)
                {
                    // Build offset distribution: ObservedApexRt - LibraryRt
                    // Use precursor lookup to get the correct library RT coordinate
                    var subsampleLookup = new Dictionary<(string, int), LibraryPrecursorInput>(subsample.Count);
                    foreach (var p in subsample)
                        if (!p.IsDecoy)
                        {
                            var key = (p.Sequence, p.ChargeState);
                            if (!subsampleLookup.ContainsKey(key))
                                subsampleLookup[key] = p;
                        }

                    var offsets = new List<double>(highQualitySub.Count);
                    foreach (var r in highQualitySub)
                    {
                        if (!subsampleLookup.TryGetValue((r.Sequence, r.ChargeState), out var p))
                            continue;
                        double libRt = p.IrtValue.HasValue ? p.IrtValue.Value
                                     : p.RetentionTime.HasValue ? p.RetentionTime.Value
                                     : double.NaN;
                        if (double.IsNaN(libRt)) continue;
                        offsets.Add(r.ObservedApexRt - libRt);
                    }

                    if (offsets.Count >= OffsetDetectionMinResults)
                    {
                        modeOffset = FindOffsetMode(offsets, OffsetKdeBandwidthMinutes);
                        offsetDetected = !double.IsNaN(modeOffset);

                        // Compute σ only on offsets near the mode to avoid
                        // outlier inflation from false matches in the broad scan
                        double coarseWindow = 1.5; // ±5 min around mode
                        var nearMode = offsets
                            .Where(o => Math.Abs(o - modeOffset) <= coarseWindow)
                            .ToList();

                        double mean = nearMode.Count > 0 ? nearMode.Average() : modeOffset;
                        double variance = nearMode.Count > 1
                            ? nearMode.Sum(o => (o - mean) * (o - mean)) / nearMode.Count
                            : 1.0;
                        double stdDev = Math.Sqrt(variance);
                        BootstrapWindowMinutes = Math.Max(3.0 * stdDev, 0.5);

                        Console.WriteLine($"  [Calibration] Offset distribution: mean={mean:F2} mode={modeOffset:F2} σ={stdDev:F2} (from {nearMode.Count}/{offsets.Count} near-mode offsets) → window=±{BootstrapWindowMinutes:F2} min");
                    }
                }

                // Step D: Full extraction centered on detected offset
                extractSw.Reset();
                extractSw.Start();

                if (offsetDetected)
                {
                    // Generate queries: LibraryRt + modeOffset ± BootstrapWindowMinutes
                    genResult = GenerateWithOffsetWindow(
                        precursors, scanIndex, baseParameters, modeOffset, BootstrapWindowMinutes);
                    Console.WriteLine($"  [Calibration] Bootstrap: offset={modeOffset:F2} min, window=±{BootstrapWindowMinutes:F1} min");
                }
                else
                {
                    // Fallback: legacy progressive widening
                    Console.WriteLine("  [Calibration] Offset detection failed, falling back to progressive widening...");
                    double[] fallbackWindows = { 1.0, 1.5, 2.0, 3.0, 5.0 };
                    bool gotEnough = false;

                    foreach (double tryWindow in fallbackWindows)
                    {
                        Console.WriteLine($"  [Calibration] Bootstrap: trying ±{tryWindow:F1} min window...");
                        var bParams = CloneWithRtTolerance(baseParameters, (float)tryWindow);
                        genResult = DiaLibraryQueryGenerator.Generate(precursors, scanIndex, bParams);
                        extractionResult = orchestrator.ExtractAll(genResult.Queries);
                        bootstrapResults = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                            precursors, genResult, extractionResult, baseParameters, scanIndex);

                        Console.WriteLine($"  [Calibration] Bootstrap ±{tryWindow:F1} min: extraction={extractSw.ElapsedMilliseconds}ms, results={bootstrapResults.Count:N0}");

                        selectSw.Start();
                        var testAnchors = SelectAnchors(bootstrapResults, precursors, genResult.PrecursorGroups,
                            InitialApexScoreThreshold, InitialTopK, requireMinFragDetRate: false);
                        selectSw.Stop();

                        Console.WriteLine($"  [Calibration] Bootstrap ±{tryWindow:F1} min: {testAnchors.Count} anchors (threshold={InitialApexScoreThreshold:F2})");

                        if (testAnchors.Count >= InitialTopK || testAnchors.Count >= MinAnchorCount * 5)
                        {
                            gotEnough = true;
                            break;
                        }
                    }

                    if (!gotEnough && (bootstrapResults == null || bootstrapResults.Count < MinAnchorCount))
                    {
                        log.Add(CreateLogEntry(0, 0, null, baseParameters.RtToleranceMinutes,
                            extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, "BootstrapFailed"));
                        return (null, null, log, bootstrapResults ?? new List<DiaSearchResult>());
                    }

                    // bootstrapResults is already set by the fallback loop
                    extractSw.Stop();

                    selectSw.Reset();
                    selectSw.Start();
                    var anchors0 = SelectAnchors(bootstrapResults, precursors, genResult.PrecursorGroups,
                        InitialApexScoreThreshold, InitialTopK, requireMinFragDetRate: false);
                    selectSw.Stop();

                    fitSw.Start();
                    currentModel = FitBootstrapModel(anchors0, out string modelLabel0);
                    fitSw.Stop();

                    if (currentModel == null)
                    {
                        log.Add(CreateLogEntry(0, anchors0.Count, null, baseParameters.RtToleranceMinutes,
                            extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, "FitFailed"));
                        return (null, null, log, bootstrapResults);
                    }

                    currentResults = bootstrapResults;
                    previousSigma = currentModel.SigmaMinutes;
                    previousSlope = currentModel.Slope;
                    double windowHW0 = Math.Max(currentModel.SigmaMinutes * SigmaMultiplier, MinWindowHalfWidthMinutes);
                    log.Add(CreateLogEntry(0, anchors0.Count, currentModel, windowHW0,
                        extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, modelLabel0));

                    // Run remaining iterations starting at 1
                    return RunRefinementIterations(
                        precursors, scanIndex, baseParameters, orchestrator,
                        currentModel, currentResults, previousSigma, previousSlope,
                        log, startIteration: 1);
                }

                // Offset-detected path: extract all precursors with offset window
                extractionResult = orchestrator.ExtractAll(genResult.Queries);
                bootstrapResults = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    precursors, genResult, extractionResult, baseParameters, scanIndex);
                extractSw.Stop();

                Console.WriteLine($"  [Calibration] Bootstrap ±{BootstrapWindowMinutes:F1} min (offset={modeOffset:F2}): extraction={extractSw.ElapsedMilliseconds}ms, results={bootstrapResults.Count:N0}");

                // Step E: Select anchors and fit initial model
                selectSw.Start();
                var anchors = SelectAnchors(bootstrapResults, precursors, genResult.PrecursorGroups,
                    InitialApexScoreThreshold, InitialTopK, requireMinFragDetRate: false);
                selectSw.Stop();

                Console.WriteLine($"  [Calibration] Bootstrap: {anchors.Count} anchors (threshold={InitialApexScoreThreshold:F2})");

                if (anchors.Count < MinAnchorCount)
                {
                    log.Add(CreateLogEntry(0, anchors.Count, null, baseParameters.RtToleranceMinutes,
                        extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, "InsufficientAnchors"));
                    return (null, null, log, bootstrapResults);
                }

                fitSw.Start();
                currentModel = FitBootstrapModel(anchors, out string modelLabel);
                fitSw.Stop();

                if (currentModel == null)
                {
                    log.Add(CreateLogEntry(0, anchors.Count, null, baseParameters.RtToleranceMinutes,
                        extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, "FitFailed"));
                    return (null, null, log, bootstrapResults);
                }

                currentResults = bootstrapResults;
                previousSigma = currentModel.SigmaMinutes;
                previousSlope = currentModel.Slope;
                double hw = Math.Max(currentModel.SigmaMinutes * SigmaMultiplier, MinWindowHalfWidthMinutes);
                log.Add(CreateLogEntry(0, anchors.Count, currentModel, hw,
                    extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, modelLabel));

                Console.WriteLine($"  [Calibration] Iteration 0: slope={currentModel.Slope:F4}, σ={currentModel.SigmaMinutes:F3}, R²={currentModel.RSquared:F4}, model={modelLabel}");
            }

            // ── Refinement Iterations 1..MaxIterations-1 ──────────────────
            return RunRefinementIterations(
                precursors, scanIndex, baseParameters, orchestrator,
                currentModel, currentResults, previousSigma, previousSlope,
                log, startIteration: 1);
        }

        // ── Refinement Loop ────────────────────────────────────────────────

        private (RtCalibrationModel Model, IRtCalibrationModel DetailedModel,
                 List<CalibrationIterationLog> Log, List<DiaSearchResult> Results)
            RunRefinementIterations(
                IList<LibraryPrecursorInput> precursors,
                DiaScanIndex scanIndex,
                DiaSearchParameters baseParameters,
                DiaExtractionOrchestrator orchestrator,
                IRtCalibrationModel currentModel,
                List<DiaSearchResult> currentResults,
                double previousSigma,
                double previousSlope,
                List<CalibrationIterationLog> log,
                int startIteration)
        {
            for (int iteration = startIteration; iteration < MaxIterations; iteration++)
            {
                Console.WriteLine($"  [Calibration] Iteration {iteration} starting...");
                var extractSw = new Stopwatch();
                var fitSw = new Stopwatch();
                var selectSw = new Stopwatch();

                extractSw.Start();

                double windowHalfWidth = Math.Max(
                    currentModel.SigmaMinutes * SigmaMultiplier,
                    MinWindowHalfWidthMinutes);

                var genResult = GenerateCalibratedWithExplicitWindow(
                    precursors, scanIndex, baseParameters, currentModel, windowHalfWidth);
                var extractionResult = orchestrator.ExtractAll(genResult.Queries);
                var results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    precursors, genResult, extractionResult, baseParameters, scanIndex);
                extractSw.Stop();

                Console.WriteLine($"  [Calibration] Iteration {iteration}: extraction={extractSw.ElapsedMilliseconds}ms, results={results.Count:N0}");

                selectSw.Start();
                var anchors = SelectAnchors(results, precursors, genResult.PrecursorGroups,
                    RefinedApexScoreThreshold, int.MaxValue, requireMinFragDetRate: true);
                selectSw.Stop();

                Console.WriteLine($"  [Calibration] Iteration {iteration}: {anchors.Count} anchors (threshold={RefinedApexScoreThreshold:F2})");

                if (anchors.Count < MinAnchorCount)
                {
                    log.Add(CreateLogEntry(iteration, anchors.Count, currentModel,
                        Math.Max(currentModel.SigmaMinutes * SigmaMultiplier, MinWindowHalfWidthMinutes),
                        extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, "InsufficientAnchors"));
                    break;
                }

                fitSw.Start();
                var libraryRts = anchors.Select(a => a.LibraryRt).ToArray();
                var observedRts = anchors.Select(a => a.ObservedRt).ToArray();

                double inlierThreshold = Math.Max(1.5 * previousSigma, 0.3);
                var linearModel = FitLinearWithRansac(libraryRts, observedRts, inlierThreshold);

                if (linearModel == null || !linearModel.IsReliable)
                {
                    log.Add(CreateLogEntry(iteration, anchors.Count, currentModel,
                        Math.Max(currentModel.SigmaMinutes * SigmaMultiplier, MinWindowHalfWidthMinutes),
                        extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, "FitFailed_Reverted"));
                    break;
                }

                IRtCalibrationModel bestModel = new LinearRtModelWrapper(linearModel);
                string modelTypeLabel = "Linear";

                if (EnableNonLinearModelSelection && anchors.Count > 200)
                {
                    bestModel = SelectBestModel(libraryRts, observedRts, linearModel, out modelTypeLabel);
                }

                fitSw.Stop();

                double newSigma = bestModel.SigmaMinutes;
                double newSlope = bestModel.Slope;

                Console.WriteLine($"  [Calibration] Iteration {iteration}: slope={newSlope:F4}, σ={newSigma:F3}, R²={bestModel.RSquared:F4}, model={modelTypeLabel}");

                // Divergence protection: more lenient in early iterations
                double divergenceThreshold = iteration <= 1 ? 1.20 : 1.10;
                if (newSigma > previousSigma * divergenceThreshold)
                {
                    Console.WriteLine($"  [Calibration] Iteration {iteration}: σ diverged ({newSigma:F3} > {previousSigma * 1.1:F3}), reverting.");
                    log.Add(CreateLogEntry(iteration, anchors.Count, bestModel,
                        Math.Max(newSigma * SigmaMultiplier, MinWindowHalfWidthMinutes),
                        extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, modelTypeLabel + "_Diverged_Reverted"));
                    break;
                }

                // Slope divergence protection
                if (iteration >= 2 && !double.IsNaN(previousSlope))
                {
                    double curDist = Math.Abs(newSlope - 1.0);
                    double prevDist = Math.Abs(previousSlope - 1.0);
                    if (curDist > prevDist + 0.02)
                    {
                        Console.WriteLine($"  [Calibration] Iteration {iteration}: slope diverging from 1.0, reverting.");
                        log.Add(CreateLogEntry(iteration, anchors.Count, bestModel,
                            Math.Max(newSigma * SigmaMultiplier, MinWindowHalfWidthMinutes),
                            extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, modelTypeLabel + "_SlopeDiverged_Reverted"));
                        break;
                    }
                }

                currentModel = bestModel;
                currentResults = results;

                double hw = Math.Max(newSigma * SigmaMultiplier, MinWindowHalfWidthMinutes);
                log.Add(CreateLogEntry(iteration, anchors.Count, bestModel, hw,
                    extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, modelTypeLabel));

                if (previousSigma > 0 && !double.IsInfinity(previousSigma))
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

        // ── Offset Detection Helpers ───────────────────────────────────────

        /// <summary>
        /// Finds the mode of an offset distribution using kernel density estimation.
        /// Evaluates KDE on a fine grid and returns the grid point with highest density.
        /// Bandwidth is fixed at OffsetKdeBandwidthMinutes.
        /// </summary>
        private static double FindOffsetMode(List<double> offsets, double bandwidth)
        {
            if (offsets.Count == 0) return double.NaN;

            double min = offsets.Min();
            double max = offsets.Max();
            if (max - min < 1e-6) return offsets[0];

            // Evaluate KDE on 500-point grid
            const int gridPoints = 500;
            double step = (max - min) / (gridPoints - 1);
            double h2 = bandwidth * bandwidth;

            double bestDensity = double.MinValue;
            double bestX = double.NaN;

            for (int g = 0; g < gridPoints; g++)
            {
                double x = min + g * step;
                double density = 0.0;
                foreach (double o in offsets)
                {
                    double d = (x - o) / bandwidth;
                    density += Math.Exp(-0.5 * d * d);
                }
                if (density > bestDensity)
                {
                    bestDensity = density;
                    bestX = x;
                }
            }

            return bestX;
        }

        /// <summary>
        /// Builds a random subsample of precursors for offset detection.
        /// Preserves the target/decoy ratio. Uses reservoir sampling.
        /// Only targets are needed for offset detection (decoys have no true RT).
        /// </summary>
        private static List<LibraryPrecursorInput> BuildSubsample(
            IList<LibraryPrecursorInput> precursors, int maxSize)
        {
            var targets = new List<LibraryPrecursorInput>(Math.Min(precursors.Count, maxSize));
            foreach (var p in precursors)
                if (!p.IsDecoy) targets.Add(p);

            if (targets.Count <= maxSize)
                return targets;

            // Fisher-Yates reservoir sample
            var rng = new Random(42);
            var reservoir = new List<LibraryPrecursorInput>(targets.Take(maxSize));
            for (int i = maxSize; i < targets.Count; i++)
            {
                int j = rng.Next(i + 1);
                if (j < maxSize)
                    reservoir[j] = targets[i];
            }

            return reservoir;
        }

        /// <summary>
        /// Generates queries with window = LibraryRt + offset ± halfWidthMinutes.
        /// Used after offset detection to center extraction on the detected shift.
        /// </summary>
        private static DiaLibraryQueryGenerator.GenerationResult GenerateWithOffsetWindow(
            IList<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            DiaSearchParameters parameters,
            double offset,
            double halfWidthMinutes)
        {
            float ppmTolerance = parameters.PpmTolerance;
            float globalRtMin = scanIndex.GetGlobalRtMin();
            float globalRtMax = scanIndex.GetGlobalRtMax();

            // First pass: count
            int totalQueryCount = 0, skippedNoWindow = 0, skippedNoFragments = 0;
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
                double libRt = p.IrtValue.HasValue ? p.IrtValue.Value
                             : p.RetentionTime.HasValue ? p.RetentionTime.Value
                             : double.NaN;

                if (!double.IsNaN(libRt))
                {
                    float center = (float)(libRt + offset);
                    // Clamp to actual run range
                    rtMin = Math.Max((float)(center - halfWidthMinutes), globalRtMin);
                    rtMax = Math.Min((float)(center + halfWidthMinutes), globalRtMax);
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

        /// <summary>
        /// Fits the bootstrap linear model from anchors, with RANSAC fallback.
        /// </summary>
        private IRtCalibrationModel FitBootstrapModel(
            List<(double LibraryRt, double ObservedRt, double QualityScore)> anchors,
            out string modelLabel)
        {
            modelLabel = "Linear";
            var libraryRts = anchors.Select(a => a.LibraryRt).ToArray();
            var observedRts = anchors.Select(a => a.ObservedRt).ToArray();

            var linearModel = FitLinearWithRansac(libraryRts, observedRts, inlierThreshold: 1.0);
            if (linearModel == null || !linearModel.IsReliable)
            {
                Console.WriteLine("  [Calibration] Iteration 0: tight RANSAC failed, retrying with 2.0 min...");
                linearModel = FitLinearWithRansac(libraryRts, observedRts, inlierThreshold: 2.0);
            }

            if (linearModel == null || !linearModel.IsReliable)
                return null;

            return new LinearRtModelWrapper(linearModel);
        }

        /// <summary>Creates a copy of parameters with a different RtToleranceMinutes.</summary>
        private static DiaSearchParameters CloneWithRtTolerance(DiaSearchParameters src, float rtTolerance)
        {
            return new DiaSearchParameters
            {
                PpmTolerance = src.PpmTolerance,
                RtToleranceMinutes = rtTolerance,
                MinFragmentsRequired = src.MinFragmentsRequired,
                MinScoreThreshold = src.MinScoreThreshold,
                MaxThreads = src.MaxThreads,
                PreferGpu = src.PreferGpu,
                CalibratedWindowSigmaMultiplier = src.CalibratedWindowSigmaMultiplier,
                ScoringStrategy = src.ScoringStrategy,
                NonlinearPower = src.NonlinearPower
            };
        }

        // ── Model Selection ────────────────────────────────────────────────

        private IRtCalibrationModel SelectBestModel(
            double[] libraryRts,
            double[] observedRts,
            RtCalibrationModel linearModel,
            out string selectedModelType)
        {
            selectedModelType = "Linear";
            IRtCalibrationModel bestModel = new LinearRtModelWrapper(linearModel);
            double bestSigma = linearModel.SigmaMinutes;

            if (linearModel.RSquared >= PiecewiseLinearRSquaredThreshold)
                return bestModel;

            try
            {
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
                // Non-linear fitting failures must never crash calibration
            }

            return bestModel;
        }

        // ── Anchor Selection ───────────────────────────────────────────────

        /// <summary>
        /// Selects high-quality calibration anchors from extraction results.
        /// Uses iRT (preferred) then RetentionTime as the library RT coordinate,
        /// matching what GenerateCalibratedWithExplicitWindow uses.
        /// </summary>
        public static List<(double LibraryRt, double ObservedRt, double QualityScore)> SelectAnchors(
            List<DiaSearchResult> results,
            IList<LibraryPrecursorInput> precursors,
            DiaLibraryQueryGenerator.PrecursorQueryGroup[] groups,
            float minApexScore,
            int maxAnchors,
            bool requireMinFragDetRate = true)
        {
            Dictionary<(string Seq, int Charge), LibraryPrecursorInput> precursorLookup = null;
            if (precursors != null && precursors.Count > 0)
            {
                precursorLookup = new Dictionary<(string, int), LibraryPrecursorInput>(precursors.Count);
                for (int i = 0; i < precursors.Count; i++)
                {
                    var p = precursors[i];
                    if (p.IsDecoy) continue;
                    var key = (p.Sequence, p.ChargeState);
                    if (!precursorLookup.ContainsKey(key))
                        precursorLookup[key] = p;
                }
            }

            var candidates = new List<(double LibraryRt, double ObservedRt, double QualityScore)>();

            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];
                if (r.IsDecoy) continue;
                if (float.IsNaN(r.ObservedApexRt)) continue;
                if (float.IsNaN(r.ApexScore) || r.ApexScore < minApexScore) continue;

                float fragDetRate = r.FragmentDetectionRate;
                if (requireMinFragDetRate && fragDetRate < 0.5f) continue;

                double libRt;
                if (precursorLookup != null &&
                    precursorLookup.TryGetValue((r.Sequence, r.ChargeState), out var precursor))
                {
                    if (precursor.IrtValue.HasValue && !double.IsNaN(precursor.IrtValue.Value))
                        libRt = precursor.IrtValue.Value;
                    else if (precursor.RetentionTime.HasValue && !double.IsNaN(precursor.RetentionTime.Value))
                        libRt = precursor.RetentionTime.Value;
                    else
                        continue;
                }
                else
                {
                    if (!r.LibraryRetentionTime.HasValue || double.IsNaN(r.LibraryRetentionTime.Value))
                        continue;
                    libRt = r.LibraryRetentionTime.Value;
                }

                double massAccuracyTerm = 1.0;
                if (!float.IsNaN(r.MeanMassErrorPpm))
                    massAccuracyTerm = 1.0 / (1.0 + Math.Abs(r.MeanMassErrorPpm) / 10.0);

                double quality = r.ApexScore * fragDetRate * massAccuracyTerm;
                candidates.Add((libRt, r.ObservedApexRt, quality));
            }

            if (candidates.Count == 0) return candidates;
            if (candidates.Count <= maxAnchors)
            {
                candidates.Sort((a, b) => b.QualityScore.CompareTo(a.QualityScore));
                return candidates;
            }

            return SelectWithRtBalance(candidates, maxAnchors, numBins: 10);
        }

        /// <summary>Legacy overload without groups.</summary>
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

        // ── Linear Model Fitting ───────────────────────────────────────────

        private static RtCalibrationModel FitLinearWithRansac(
            double[] libraryRts, double[] observedRts, double inlierThreshold)
        {
            var fitOptions = new RtCalibrationFitter.FitOptions
            {
                UseRansac = true,
                RansacInlierThresholdMinutes = inlierThreshold,
                RansacMinInlierFraction = 0.5
            };

            return RtCalibrationFitter.Fit(
                (IList<double>)libraryRts,
                (IList<double>)observedRts,
                fitOptions);
        }

        // ── Calibrated Query Generation ────────────────────────────────────

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

            int totalQueryCount = 0, skippedNoWindow = 0, skippedNoFragments = 0;
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

            int queryIndex = 0;
            for (int i = 0; i < precursors.Count; i++)
            {
                var p = precursors[i];
                int windowId = scanIndex.FindWindowForPrecursorMz(p.PrecursorMz);
                if (windowId < 0 || p.FragmentCount == 0)
                    continue;

                float rtMin, rtMax;

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