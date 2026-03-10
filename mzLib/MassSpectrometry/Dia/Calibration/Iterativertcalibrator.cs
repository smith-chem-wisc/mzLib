// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// Phase 15, Prompts 2-3 (fixed in Prompt 5): Iterative calibration framework with model selection
// Phase 24: T/D ratio scan bootstrap — works for any library scale (iRT or native RT).
// Placement: MassSpectrometry/Dia/Calibration/IterativeRtCalibrator.cs

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using MassSpectrometry.Dia.Calibration;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Implements iterative RT calibration with a T/D ratio scan bootstrap.
    ///
    /// Bootstrap (Iteration 0):
    ///   1. Extract all precursors over the full run RT range in one pass.
    ///   2. Divide library RT into 1-minute bins.
    ///   3. Within each bin, compute residual = observedApexRt − libRt, then
    ///      bucket residuals into 20-second increments and find the bucket with
    ///      the highest T/D ratio.  That bucket is the true offset for that
    ///      slice of library RT space.
    ///   4. Fit a line through all (libRt_mid, libRt_mid + peak_offset) control
    ///      points → bootstrap calibration model.
    ///
    /// This works for any library type:
    ///   - Native-RT libraries (slope ≈ 1, offset ≈ 0): all bins agree on ~0 offset.
    ///   - iRT-scaled libraries (slope ≈ 0.18): each bin gives a different offset,
    ///     tracing the correct slope across the library RT range.
    ///
    /// Refinement (Iterations 1+):
    ///   - Use the bootstrap model to generate per-precursor narrow RT windows.
    ///   - Select high-quality anchors, refit model, iterate until convergence.
    ///
    /// Model Selection (refinement only):
    ///   - Iterations 0–2: always Linear.
    ///   - Iteration 3+: if linear R² &lt; 0.995 and anchors &gt; 200, try PiecewiseLinear.
    ///   - LOWESS tried only if PiecewiseLinear R² &lt; 0.990.
    /// </summary>
    ///    a. Extract a random subsample of precursors using the full run RT range
    ///    b. Score results, compute per-result offset = ObservedApexRt - LibraryRt
    public class IterativeRtCalibrator
    {
        // ── Configuration ──────────────────────────────────────────────────

        /// <summary>Maximum iterations including bootstrap. Default 4.</summary>
        public int MaxIterations { get; set; } = 6;

        /// <summary>Stop when Δσ/σ &lt; this threshold (5%). Default 0.05.</summary>
        public double ConvergenceThreshold { get; set; } = 0.02;

        /// <summary>
        /// Optional progress reporter for calibration log messages.
        /// When set, messages are routed here instead of Console.WriteLine.
        /// Used by MetaMorpheus to surface calibration progress in the GUI/log.
        /// </summary>
        public Action<string> ProgressReporter { get; set; }

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

            // ── Bootstrap: Center-Outward T/D Ratio Scan ─────────────────
            // Start with the middle 20% of library RT (iRT ranks 40–60%).  Those
            // precursors elute near the gradient midpoint regardless of slope, so a
            // flat RT window centred on the TIC peak captures them cleanly.  Run the
            // T/D scan on just that band to get an initial slope estimate, optionally
            // tighten it with one refinement pass, then walk outward band by band
            // (20–40 / 60–80, then 0–20 / 80–100), using the current model to narrow
            // the window for each successive band.
            {
                var bootstrapSw = Stopwatch.StartNew();

                // ── Sort all targets by library RT ───────────────────────────
                float globalRtMin = scanIndex.GetGlobalRtMin();
                float globalRtMax = scanIndex.GetGlobalRtMax();

                var targets = new List<(int PrecIdx, double LibRt)>(precursors.Count);
                for (int i = 0; i < precursors.Count; i++)
                {
                    var p = precursors[i];
                    if (p.IsDecoy) continue;
                    double lr = p.IrtValue.HasValue ? p.IrtValue.Value
                              : p.RetentionTime.HasValue ? p.RetentionTime.Value
                              : double.NaN;
                    if (!double.IsNaN(lr)) targets.Add((i, lr));
                }
                targets.Sort((a, b) => a.LibRt.CompareTo(b.LibRt));

                if (targets.Count < MinAnchorCount * 2)
                {
                    Console.WriteLine($"  [Calibration] Too few targets ({targets.Count}), skipping T/D scan bootstrap.");
                    return LegacyProgressiveBootstrap(precursors, scanIndex, baseParameters,
                        orchestrator, log, ref currentModel, ref currentResults,
                        ref previousSigma, ref previousSlope);
                }

                int n = targets.Count;
                // Five equal bands by sorted iRT rank: 0–20, 20–40, 40–60, 60–80, 80–100%
                int[] bandLo = { 0, n / 5, 2 * n / 5, 3 * n / 5, 4 * n / 5 };
                int[] bandHi = { n / 5, 2 * n / 5, 3 * n / 5, 4 * n / 5, n };

                // ── TIC midpoint for center-band flat window ──────────────────
                var (elutionMin, elutionMax) = EstimatePeptideElutionWindow(scanIndex);
                double elutionMid = (elutionMin + elutionMax) / 2.0;
                // Initial flat half-width: 20% of run length, min 3 min, max 8 min
                double centerHalfWin = Math.Clamp((elutionMax - elutionMin) * 0.20, 3.0, 8.0);

                Console.WriteLine($"  [Calibration] Gradient window: {elutionMin:F1}–{elutionMax:F1} min (mid={elutionMid:F1})");
                Console.WriteLine($"  [Calibration] Bootstrap: {n} targets, center band = iRT ranks 40–60%");

                // Accumulated control points across all bands
                var allControlPoints = new List<(double LibRtMid, double BestOffset)>(32);

                // ── Band 2 (center, ranks 40–60%): flat RT window ─────────────
                {
                    var bandIndices = GetBandPrecursorIndices(targets, bandLo[2], bandHi[2]);
                    var bandPrecursors = BuildBandPrecursorList(precursors, bandIndices, includeDecoys: true);

                    var sw = Stopwatch.StartNew();
                    var gen = GenerateWithFlatWindowForAll(bandPrecursors, scanIndex, baseParameters,
                        (float)(elutionMid - centerHalfWin), (float)(elutionMid + centerHalfWin));
                    var ext = orchestrator.ExtractAll(gen.Queries);
                    var res = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                        bandPrecursors, gen, ext, baseParameters, scanIndex);
                    sw.Stop();

                    var pts = TdRatioScanBand(res, bandPrecursors, gen);
                    Console.WriteLine($"  [Calibration] Center band: {bandPrecursors.Count} precursors, " +
                        $"{pts.Count} control points ({sw.ElapsedMilliseconds}ms)");
                    allControlPoints.AddRange(pts);

                    // One tightening pass if we got a decent model
                    var modelCenter = FitFromControlPoints(allControlPoints);
                    if (modelCenter != null && pts.Count >= 3)
                    {
                        double tightWin = Math.Max(modelCenter.SigmaMinutes * 3.0, MinWindowHalfWidthMinutes);
                        sw.Restart();
                        gen = GenerateCalibratedWithExplicitWindow(bandPrecursors, scanIndex,
                            baseParameters, modelCenter, tightWin);
                        ext = orchestrator.ExtractAll(gen.Queries);
                        res = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                            bandPrecursors, gen, ext, baseParameters, scanIndex);
                        sw.Stop();

                        var pts2 = TdRatioScanBand(res, bandPrecursors, gen);
                        Console.WriteLine($"  [Calibration] Center band (tightened ±{tightWin:F2} min): " +
                            $"{pts2.Count} control points ({sw.ElapsedMilliseconds}ms)");
                        // Replace center-band points with the tighter ones if better
                        if (pts2.Count >= pts.Count)
                        {
                            allControlPoints.Clear();
                            allControlPoints.AddRange(pts2);
                        }
                    }
                }

                // ── Bands 1+3 (mid, ranks 20–40% and 60–80%) ─────────────────
                {
                    var modelSoFar = FitFromControlPoints(allControlPoints);
                    if (modelSoFar != null)
                    {
                        // σ×4 from the tightened center model, but floor at 1.5 min to
                        // absorb extrapolation error from projecting center-band model outward.
                        double winB = Math.Max(modelSoFar.SigmaMinutes * 4.0, 1.5);
                        var bandIndices = GetBandPrecursorIndices(targets,
                            new[] { (bandLo[1], bandHi[1]), (bandLo[3], bandHi[3]) });
                        var bandPrecursors = BuildBandPrecursorList(precursors, bandIndices, includeDecoys: true);

                        var sw = Stopwatch.StartNew();
                        var gen = GenerateCalibratedWithExplicitWindow(bandPrecursors, scanIndex,
                            baseParameters, modelSoFar, winB);
                        var ext = orchestrator.ExtractAll(gen.Queries);
                        var res = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                            bandPrecursors, gen, ext, baseParameters, scanIndex);
                        sw.Stop();

                        var pts = TdRatioScanBand(res, bandPrecursors, gen);
                        Console.WriteLine($"  [Calibration] Mid bands (±{winB:F2} min): " +
                            $"{bandPrecursors.Count} precursors, {pts.Count} control points ({sw.ElapsedMilliseconds}ms)");
                        allControlPoints.AddRange(pts);
                    }
                    else
                    {
                        Console.WriteLine("  [Calibration] Center band yielded no model — skipping mid bands.");
                    }
                }

                // ── Bands 0+4 (outer, ranks 0–20% and 80–100%) ───────────────
                {
                    var modelSoFar = FitFromControlPoints(allControlPoints);
                    if (modelSoFar != null)
                    {
                        double winC = Math.Max(modelSoFar.SigmaMinutes * 4.0, 2.0);
                        var bandIndices = GetBandPrecursorIndices(targets,
                            new[] { (bandLo[0], bandHi[0]), (bandLo[4], bandHi[4]) });
                        var bandPrecursors = BuildBandPrecursorList(precursors, bandIndices, includeDecoys: true);

                        var sw = Stopwatch.StartNew();
                        var gen = GenerateCalibratedWithExplicitWindow(bandPrecursors, scanIndex,
                            baseParameters, modelSoFar, winC);
                        var ext = orchestrator.ExtractAll(gen.Queries);
                        var res = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                            bandPrecursors, gen, ext, baseParameters, scanIndex);
                        sw.Stop();

                        var pts = TdRatioScanBand(res, bandPrecursors, gen);
                        Console.WriteLine($"  [Calibration] Outer bands (±{winC:F2} min): " +
                            $"{bandPrecursors.Count} precursors, {pts.Count} control points ({sw.ElapsedMilliseconds}ms)");
                        allControlPoints.AddRange(pts);
                    }
                }

                // ── Fit final bootstrap model from all control points ─────────
                // First remove control points whose offset deviates more than
                // 3× MAD from the median offset — catches the tail-of-gradient
                // early-eluters that inflate σ and widen the full-pass window.
                var filteredControlPoints = FilterOutlierControlPoints(allControlPoints);
                Console.WriteLine($"  [Calibration] Bootstrap control points: {allControlPoints.Count} total, " +
                    $"{filteredControlPoints.Count} after outlier filter");

                var bootstrapModel = FitFromControlPoints(filteredControlPoints)
                               ?? FitFromControlPoints(allControlPoints); // fallback to unfiltered
                bootstrapSw.Stop();

                if (bootstrapModel == null)
                {
                    Console.WriteLine("  [Calibration] Bootstrap T/D scan produced no model — falling back.");
                    log.Add(CreateLogEntry(0, 0, null, baseParameters.RtToleranceMinutes,
                        bootstrapSw.Elapsed, TimeSpan.Zero, TimeSpan.Zero, "TdScan_NoModel_Fallback"));
                    return LegacyProgressiveBootstrap(precursors, scanIndex, baseParameters,
                        orchestrator, log, ref currentModel, ref currentResults,
                        ref previousSigma, ref previousSlope);
                }

                Console.WriteLine($"  [Calibration] Bootstrap complete: slope={bootstrapModel.Slope:F4}, " +
                    $"σ={bootstrapModel.SigmaMinutes:F3}, R²={bootstrapModel.RSquared:F4} " +
                    $"({filteredControlPoints.Count} control points, {bootstrapSw.ElapsedMilliseconds}ms)");

                // ── Full extraction with the bootstrap model ──────────────────
                double fullWin = Math.Max(bootstrapModel.SigmaMinutes * SigmaMultiplier, MinWindowHalfWidthMinutes);
                var swFull = Stopwatch.StartNew();
                var genFull = GenerateCalibratedWithExplicitWindow(
                    precursors, scanIndex, baseParameters, bootstrapModel, fullWin);
                var extFull = orchestrator.ExtractAll(genFull.Queries);
                currentResults = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    precursors, genFull, extFull, baseParameters, scanIndex);
                swFull.Stop();

                Console.WriteLine($"  [Calibration] Bootstrap full pass: {currentResults.Count:N0} results, window ±{fullWin:F2} min ({swFull.ElapsedMilliseconds}ms)");

                log.Add(CreateLogEntry(0, filteredControlPoints.Count, bootstrapModel, fullWin,
                    bootstrapSw.Elapsed, TimeSpan.Zero, swFull.Elapsed, "TdScan_Linear"));
                currentModel = bootstrapModel;
                previousSigma = bootstrapModel.SigmaMinutes;
                previousSlope = bootstrapModel.Slope;
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
        /// <summary>
        /// Generates queries with a fixed [rtMin, rtMax] window for every precursor,
        /// regardless of library RT. Used in the T/D scan bootstrap.
        /// </summary>
        private static DiaLibraryQueryGenerator.GenerationResult GenerateWithFlatWindowForAll(
            IList<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            DiaSearchParameters parameters,
            float rtMin,
            float rtMax)
        {
            float ppmTolerance = parameters.PpmTolerance;
            rtMin = Math.Max(rtMin, scanIndex.GetGlobalRtMin());
            rtMax = Math.Min(rtMax, scanIndex.GetGlobalRtMax());

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
                if (windowId < 0 || p.FragmentCount == 0) continue;

                int groupOffset = queryIndex;
                for (int f = 0; f < p.FragmentCount; f++)
                {
                    queries[queryIndex++] = new FragmentQuery(
                        targetMz: p.FragmentMzs[f],
                        tolerancePpm: ppmTolerance,
                        rtMin: rtMin,
                        rtMax: rtMax,
                        windowId: windowId,
                        queryId: queryIndex - 1);
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
        /// Estimates the peptide elution window using TIC (sum of scan intensities).
        /// Smooths TIC with a 5-scan median filter, finds the contiguous span where
        /// smoothed TIC ≥ 10% of peak TIC.  Falls back to full run bounds on failure.
        /// </summary>
        private (double ElutionMin, double ElutionMax) EstimatePeptideElutionWindow(DiaScanIndex scanIndex)
        {
            int nScans = scanIndex.ScanCount;
            if (nScans < 5) return (scanIndex.GetGlobalRtMin(), scanIndex.GetGlobalRtMax());

            var tics = new double[nScans];
            for (int s = 0; s < nScans; s++)
            {
                var intensities = scanIndex.GetScanIntensitySpan(s);
                double sum = 0;
                foreach (float v in intensities) sum += v;
                tics[s] = sum;
            }

            // 5-scan median smooth
            var smooth = new double[nScans];
            for (int s = 0; s < nScans; s++)
            {
                int lo = Math.Max(0, s - 2), hi = Math.Min(nScans - 1, s + 2);
                int len = hi - lo + 1;
                var buf = new double[len];
                for (int k = 0; k < len; k++) buf[k] = tics[lo + k];
                Array.Sort(buf);
                smooth[s] = buf[len / 2];
            }

            double peakTic = 0;
            foreach (double v in smooth) if (v > peakTic) peakTic = v;
            if (peakTic <= 0) return (scanIndex.GetGlobalRtMin(), scanIndex.GetGlobalRtMax());

            double threshold = peakTic * 0.10;
            double minRt = double.MaxValue, maxRt = double.MinValue;
            for (int s = 0; s < nScans; s++)
            {
                if (smooth[s] >= threshold)
                {
                    double rt = scanIndex.GetScanRt(s);
                    if (rt < minRt) minRt = rt;
                    if (rt > maxRt) maxRt = rt;
                }
            }

            if (minRt >= maxRt)
                return (scanIndex.GetGlobalRtMin(), scanIndex.GetGlobalRtMax());

            return (minRt, maxRt);
        }

        /// <summary>
        /// Builds a subset of precursors for a bootstrap band.
        /// Always includes all decoys (for FDR scoring), plus the targets
        /// whose input indices are listed in targetIndices.
        /// </summary>
        private static List<LibraryPrecursorInput> BuildBandPrecursorList(
            IList<LibraryPrecursorInput> all,
            IList<int> targetIndices,
            bool includeDecoys)
        {
            var result = new List<LibraryPrecursorInput>(targetIndices.Count);
            foreach (int idx in targetIndices) result.Add(all[idx]);
            if (includeDecoys)
                for (int i = 0; i < all.Count; i++)
                    if (all[i].IsDecoy) result.Add(all[i]);
            return result;
        }

        /// <summary>
        /// Legacy progressive-widening bootstrap fallback.
        /// Tries ±1, ±1.5, ±2, ±3, ±5 min until enough anchors are found.
        /// </summary>
        private (RtCalibrationModel Model, IRtCalibrationModel DetailedModel,
                 List<CalibrationIterationLog> Log, List<DiaSearchResult> Results)
            LegacyProgressiveBootstrap(
                IList<LibraryPrecursorInput> precursors,
                DiaScanIndex scanIndex,
                DiaSearchParameters baseParameters,
                DiaExtractionOrchestrator orchestrator,
                List<CalibrationIterationLog> log,
                ref IRtCalibrationModel currentModel,
                ref List<DiaSearchResult> currentResults,
                ref double previousSigma,
                ref double previousSlope)
        {
            Console.WriteLine("  [Calibration] Legacy progressive-widening bootstrap...");
            double[] fallbackWindows = { 1.0, 1.5, 2.0, 3.0, 5.0 };
            DiaLibraryQueryGenerator.GenerationResult genResult = default;
            List<DiaSearchResult> bootstrapResults = null;

            var extractSw = new Stopwatch();
            var fitSw = new Stopwatch();
            var selectSw = new Stopwatch();

            foreach (double tryWindow in fallbackWindows)
            {
                Console.WriteLine($"  [Calibration] Bootstrap: trying ±{tryWindow:F1} min window...");
                extractSw.Restart();
                var bParams = CloneWithRtTolerance(baseParameters, (float)tryWindow);
                genResult = DiaLibraryQueryGenerator.Generate(precursors, scanIndex, bParams);
                var extractionResult = orchestrator.ExtractAll(genResult.Queries);
                bootstrapResults = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    precursors, genResult, extractionResult, baseParameters, scanIndex);
                extractSw.Stop();

                Console.WriteLine($"  [Calibration] Bootstrap ±{tryWindow:F1} min: {bootstrapResults.Count:N0} results ({extractSw.ElapsedMilliseconds}ms)");

                selectSw.Restart();
                var testAnchors = SelectAnchors(bootstrapResults, precursors, genResult.PrecursorGroups,
                    InitialApexScoreThreshold, InitialTopK, requireMinFragDetRate: false);
                selectSw.Stop();

                Console.WriteLine($"  [Calibration] Bootstrap ±{tryWindow:F1} min: {testAnchors.Count} anchors");

                if (testAnchors.Count >= MinAnchorCount)
                {
                    fitSw.Restart();
                    currentModel = FitBootstrapModel(testAnchors, out string ml);
                    fitSw.Stop();

                    if (currentModel != null)
                    {
                        currentResults = bootstrapResults;
                        previousSigma = currentModel.SigmaMinutes;
                        previousSlope = currentModel.Slope;
                        double hw = Math.Max(currentModel.SigmaMinutes * SigmaMultiplier, MinWindowHalfWidthMinutes);
                        log.Add(CreateLogEntry(0, testAnchors.Count, currentModel, hw,
                            extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, "Legacy_" + ml));
                        Console.WriteLine($"  [Calibration] Legacy bootstrap: slope={currentModel.Slope:F4}, σ={currentModel.SigmaMinutes:F3}");
                        break;
                    }
                }
            }

            if (currentModel == null)
            {
                log.Add(CreateLogEntry(0, 0, null, baseParameters.RtToleranceMinutes,
                    extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, "BootstrapFailed"));
                return (null, null, log, bootstrapResults ?? new List<DiaSearchResult>());
            }

            return RunRefinementIterations(
                precursors, scanIndex, baseParameters, orchestrator,
                currentModel, currentResults, previousSigma, previousSlope,
                log, startIteration: 1);
        }

        /// <summary>Returns precursor indices for a single contiguous iRT rank band [lo, hi).</summary>
        private static List<int> GetBandPrecursorIndices(
            List<(int PrecIdx, double LibRt)> targets, int lo, int hi)
        {
            var idx = new List<int>(hi - lo);
            for (int i = lo; i < hi; i++) idx.Add(targets[i].PrecIdx);
            return idx;
        }

        /// <summary>Returns precursor indices for multiple (lo, hi) band ranges.</summary>
        private static List<int> GetBandPrecursorIndices(
            List<(int PrecIdx, double LibRt)> targets,
            IEnumerable<(int Lo, int Hi)> bands)
        {
            var idx = new List<int>();
            foreach (var (lo, hi) in bands)
                for (int i = lo; i < hi; i++) idx.Add(targets[i].PrecIdx);
            return idx;
        }

        /// <summary>
        /// T/D ratio scan over one band's extraction results.
        /// Divides library RT into 1-minute slices; within each slice bins the residual
        /// (observedApexRt − libRt) into 20-second buckets and finds the bucket with
        /// the highest target/decoy ratio.  Returns one (libRtMid, peakOffset) control
        /// point per slice that passes the quality threshold (≥ 3 results, T/D ≥ 2:1).
        /// </summary>
        private List<(double LibRtMid, double BestOffset)> TdRatioScanBand(
            List<DiaSearchResult> results,
            IList<LibraryPrecursorInput> bandPrecursors,
            DiaLibraryQueryGenerator.GenerationResult gen)
        {
            const double libBinWidth = 1.0;
            const double residualBinMin = 20.0 / 60.0;  // 20-second buckets
            const int minBucketTotal = 3;
            const double minTdRatio = 2.0;

            var data = new List<(double LibRt, double ObsRt, bool IsDecoy)>(results.Count);
            for (int gi = 0; gi < gen.PrecursorGroups.Length && gi < results.Count; gi++)
            {
                var r = results[gi];
                if (float.IsNaN(r.ObservedApexRt) || r.ObservedApexRt <= 0f) continue;
                var prec = bandPrecursors[gen.PrecursorGroups[gi].InputIndex];
                double lr = prec.IrtValue.HasValue ? prec.IrtValue.Value
                          : prec.RetentionTime.HasValue ? prec.RetentionTime.Value
                          : double.NaN;
                if (!double.IsNaN(lr))
                    data.Add((lr, r.ObservedApexRt, prec.IsDecoy));
            }

            if (data.Count < MinAnchorCount) return new List<(double, double)>();

            double libMin = double.MaxValue, libMax = double.MinValue;
            double resMin = double.MaxValue, resMax = double.MinValue;
            foreach (var d in data)
            {
                if (d.LibRt < libMin) libMin = d.LibRt;
                if (d.LibRt > libMax) libMax = d.LibRt;
                double res = d.ObsRt - d.LibRt;
                if (res < resMin) resMin = res;
                if (res > resMax) resMax = res;
            }
            if (resMax - resMin < residualBinMin) return new List<(double, double)>();

            int numLibBins = Math.Max(1, (int)Math.Ceiling((libMax - libMin) / libBinWidth));
            int numResBins = (int)Math.Ceiling((resMax - resMin) / residualBinMin) + 1;

            var tgt = new int[numLibBins, numResBins];
            var dec = new int[numLibBins, numResBins];
            foreach (var d in data)
            {
                int lb = Math.Min((int)((d.LibRt - libMin) / libBinWidth), numLibBins - 1);
                int rb = Math.Clamp((int)((d.ObsRt - d.LibRt - resMin) / residualBinMin), 0, numResBins - 1);
                if (d.IsDecoy) dec[lb, rb]++; else tgt[lb, rb]++;
            }

            var pts = new List<(double LibRtMid, double BestOffset)>(numLibBins);
            for (int lb = 0; lb < numLibBins; lb++)
            {
                double bestRatio = -1.0; int bestRb = -1;
                for (int rb = 0; rb < numResBins; rb++)
                {
                    if (tgt[lb, rb] + dec[lb, rb] < minBucketTotal) continue;
                    double ratio = (double)tgt[lb, rb] / (dec[lb, rb] + 1.0);
                    if (ratio > bestRatio) { bestRatio = ratio; bestRb = rb; }
                }
                if (bestRb < 0 || tgt[lb, bestRb] < 2 || bestRatio < minTdRatio) continue;

                double offset = resMin + (bestRb + 0.5) * residualBinMin;
                double libMid = libMin + (lb + 0.5) * libBinWidth;
                pts.Add((libMid, offset));
                Console.WriteLine($"    libRt={libMid:F1}  offset={offset:+0.000;-0.000} min  " +
                    $"obsRt≈{libMid + offset:F2}  T/D={tgt[lb, bestRb]}/{dec[lb, bestRb]}");
            }
            return pts;
        }

        /// <summary>
        /// Fits a linear calibration model from accumulated (libRtMid, peakOffset)
        /// control points via OLS+RANSAC.  Each point contributes
        /// (libRtMid, libRtMid + peakOffset) as an (x, y) anchor.
        /// Returns null if fewer than 3 points or the fit is unreliable.
        /// </summary>
        /// <summary>
        /// Removes control points whose offset deviates more than 3× MAD from the
        /// median offset.  This eliminates late-gradient early-eluters and other
        /// systematic outlier bins that inflate σ and widen the full-pass window.
        /// Always retains at least 3 points; falls back to all points if fewer remain.
        /// </summary>
        private static List<(double LibRtMid, double BestOffset)> FilterOutlierControlPoints(
            List<(double LibRtMid, double BestOffset)> pts)
        {
            if (pts.Count <= 3) return pts;

            // Compute median offset
            var offsets = pts.Select(p => p.BestOffset).OrderBy(x => x).ToArray();
            double median = offsets[offsets.Length / 2];

            // MAD = median of |offset - median|
            var absDevs = offsets.Select(o => Math.Abs(o - median)).OrderBy(x => x).ToArray();
            double mad = absDevs[absDevs.Length / 2];

            // Use 1.4826 × MAD as robust σ estimate; threshold at 3×
            double robustSigma = 1.4826 * mad;
            double threshold = Math.Max(robustSigma * 3.0, 0.5); // floor at 0.5 min

            var filtered = pts.Where(p => Math.Abs(p.BestOffset - median) <= threshold).ToList();
            return filtered.Count >= 3 ? filtered : pts;
        }

        private IRtCalibrationModel FitFromControlPoints(
            List<(double LibRtMid, double BestOffset)> pts)
        {
            if (pts.Count < 3) return null;

            var xs = new double[pts.Count];
            var ys = new double[pts.Count];
            for (int i = 0; i < pts.Count; i++)
            {
                xs[i] = pts[i].LibRtMid;
                ys[i] = pts[i].LibRtMid + pts[i].BestOffset;
            }

            // Control points are already T/D-verified — for small sets (≤10) the
            // x-range may be narrow (center band only spans ~5 min) making R²
            // unreliable, so use direct OLS.  RANSAC for larger sets.
            RtCalibrationFitter.FitOptions fitOpts;
            if (pts.Count <= 10)
            {
                fitOpts = new RtCalibrationFitter.FitOptions
                {
                    UseRansac = false,
                    OutlierRejectionPasses = 1,
                    MinAnchors = 3
                };
            }
            else
            {
                fitOpts = new RtCalibrationFitter.FitOptions
                {
                    UseRansac = true,
                    RansacInlierThresholdMinutes = 1.5,
                    RansacMinInlierFraction = 0.5,
                    OutlierRejectionPasses = 2,
                    MinAnchors = 3
                };
            }

            var model = RtCalibrationFitter.Fit(
                (ReadOnlySpan<double>)xs, (ReadOnlySpan<double>)ys, fitOpts);

            // Accept any finite-slope model — IsReliable can reject good fits when
            // x-range is narrow (center band ~5 min wide gives artificially low R²).
            return model == null ? null : new LinearRtModelWrapper(model);
        }

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
                (ReadOnlySpan<double>)libraryRts,
                (ReadOnlySpan<double>)observedRts,
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