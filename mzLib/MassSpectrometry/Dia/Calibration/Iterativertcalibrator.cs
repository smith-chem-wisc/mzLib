// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// Phase 15, Prompts 2-3 (fixed in Prompt 5): Iterative calibration framework with model selection
// Phase 21: Offset-detection bootstrap (replaced — only worked when slope ≈ 1).
// Phase 23: Center-outward bootstrap — works for any iRT scale (slope ≠ 1).
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
    /// Algorithm (Phase 23 — center-outward bootstrap):
    /// 1. Estimate the peptide elution window from the MS2 TIC profile.
    /// 2. Sort target precursors by iRT rank (not value). Extract the center quintile
    ///    (rank 40–60%) using a flat RT window around the TIC midpoint. Fit model_0.
    ///    This works regardless of iRT scale because center-iRT peptides always elute
    ///    near the gradient midpoint.
    /// 3. Expand outward: extract mid bands (rank 20–40%, 60–80%) using model_0 ± 2σ.
    ///    Fit model_1 from all anchors so far.
    /// 4. Expand to outer bands (rank 0–20%, 80–100%) using model_1 ± 2σ. Fit model_2.
    /// 5. Full bootstrap pass: all precursors using model_2 ± 4σ.
    /// 6. Iterative refinement until convergence.
    ///
    /// This replaces the Phase 21 offset-detection bootstrap, which assumed a single
    /// near-constant iRT→RT offset (slope ≈ 1). That assumption breaks for Koina/Prosit
    /// libraries where the true slope is ~0.175 and the offset varies by 167 min across
    /// the iRT range, making KDE mode detection meaningless.
    ///
    /// Model Selection Logic:
    /// - Phases A/B/C: always Linear (building the model incrementally)
    /// - Refinement iteration 1+: if linear R² &lt; 0.995 and anchors &gt; 200, try PiecewiseLinear
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

        /// <summary>Window = predicted ± N × σ. Default 2.0 — with true σ ≈ 0.15–0.25 min
        /// this gives a ±0.3–0.5 min window, matching the expected extraction width.
        /// Note: the calibration window cap uses BootstrapSigmaMultiplier × bootstrapSigma,
        /// so this value only affects the final search window via DiaCalibrationPipeline.</summary>
        public double SigmaMultiplier { get; set; } = 2.0;

        /// <summary>Maximum anchors for bootstrap phases (A/B/C). Uncapped — more anchors = better slope estimate.</summary>
        public int InitialTopK { get; set; } = int.MaxValue;

        /// <summary>Strict ApexScore threshold for bootstrap iteration. Default 0.85.</summary>
        public float InitialApexScoreThreshold { get; set; } = 0.85f;

        /// <summary>Relaxed ApexScore threshold for later iterations. Default 0.5.</summary>
        public float RefinedApexScoreThreshold { get; set; } = 0.5f;

        /// <summary>Absolute floor on any window half-width. Default 0.3 min.</summary>
        public double MinWindowHalfWidthMinutes { get; set; } = 0.3;

        /// <summary>Floor on refinement extraction window. Default 0.5 min — prevents
        /// window from going below one scan cycle's worth of RT coverage.</summary>
        public double MinRefinementWindowMinutes { get; set; } = 0.5;

        /// <summary>Minimum anchors required to proceed with fitting. Default 20.</summary>
        public int MinAnchorCount { get; set; } = 20;

        /// <summary>Precursors to subsample for offset detection. Default 5000.</summary>
        public int OffsetDetectionSubsampleSize { get; set; } = 5000;

        /// <summary>Half-width of extraction window after offset detection. Default 2.0 min.</summary>
        public double BootstrapWindowMinutes { get; set; } = 2.0;

        /// <summary>Minimum high-quality subsample results to trust the detected mode. Default 50.</summary>
        public int OffsetDetectionMinResults { get; set; } = 50;

        /// <summary>KDE bandwidth for offset mode detection. Default 0.5 min.</summary>
        public double OffsetKdeBandwidthMinutes { get; set; } = 0.5;

        /// <summary>
        /// Half-width of the flat RT window used for Phase A (center-band) extraction.
        /// Centered on the TIC midpoint. Wide enough to catch center-iRT peptides
        /// regardless of unknown slope. Default 4.0 min.
        /// </summary>
        public double CenterBandHalfWidthMinutes { get; set; } = 4.0;

        /// <summary>
        /// Sigma multiplier used for Phase B and C bootstrap windows.
        /// Tighter than the main refinement multiplier to keep the expanding
        /// windows from growing too wide before the model is well-constrained.
        /// Default 2.0.
        /// </summary>
        public double BootstrapSigmaMultiplier { get; set; } = 2.0;

        /// <summary>
        /// Fraction of peak TIC used as the threshold for detecting the peptide
        /// elution window from the MS2 TIC profile. Default 0.10 (10%).
        /// </summary>
        public double TicThresholdFraction { get; set; } = 0.10;

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

        // ── Progress Reporting ─────────────────────────────────────────────

        /// <summary>
        /// Optional callback invoked for every calibration status message.
        /// When set (e.g. to MetaMorpheusEngine.Status), messages are surfaced
        /// to the MetaMorpheus GUI/log in addition to Console.WriteLine.
        /// </summary>
        public Action<string> ProgressReporter { get; set; }

        /// <summary>
        /// Writes a calibration message to Console and Debug output, and if set,
        /// also forwards it to ProgressReporter so MetaMorpheus can display it.
        /// </summary>
        private void Report(string message)
        {
            Console.WriteLine(message);
            System.Diagnostics.Debug.WriteLine(message);
            ProgressReporter?.Invoke(message);
        }

        /// <summary>
        /// Reports the target/decoy ratio and counts from a result list.
        /// Also reports how many pass the apex score threshold.
        /// Called after every extraction so we can see T/D improving across phases.
        /// </summary>
        private void ReportTdRatio(string label, List<DiaSearchResult> results, float apexThreshold)
        {
            int targets = 0, decoys = 0, targetsAbove = 0, decoysAbove = 0;
            foreach (var r in results)
            {
                bool above = !float.IsNaN(r.ApexScore) && r.ApexScore >= apexThreshold;
                if (r.IsDecoy) { decoys++; if (above) decoysAbove++; }
                else { targets++; if (above) targetsAbove++; }
            }
            double allRatio = decoys > 0 ? (double)targets / decoys : double.PositiveInfinity;
            double aboveRatio = decoysAbove > 0 ? (double)targetsAbove / decoysAbove : double.PositiveInfinity;
            Report($"  [Calibration] {label}: {targets}T / {decoys}D (ratio {allRatio:F2}) | " +
                   $"ApexScore≥{apexThreshold:F2}: {targetsAbove}T / {decoysAbove}D (ratio {aboveRatio:F2})");
        }

        // ── Main Entry Point ───────────────────────────────────────────────

        /// <summary>
        /// Runs the full iterative calibration loop.
        ///
        /// Bootstrap (iteration 0): offset-detection.
        ///   a. Extract a subsample over the full RT range.
        ///   b. For high-quality results, compute offset = ObservedApexRt - LibraryRt.
        ///   c. Find the mode of the offset distribution via KDE — targets cluster at
        ///      the true shift; decoy offsets are flat/random so they don't contribute
        ///      a mode. This is immune to the T/D=1 raw-count problem.
        ///   d. Re-extract ALL precursors at LibraryRt + modeOffset ± window.
        ///   e. Fit initial linear model.
        /// Iterations 1-N: re-extract at predicted RT ± Nσ, refit, until convergence.
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

            var extractSw = new Stopwatch();
            var fitSw = new Stopwatch();
            var selectSw = new Stopwatch();

            Report("  [Calibration] Iteration 0: center-outward bootstrap (Phase 23)...");

            // ── Estimate elution window from file RT range ─────────────────────────
            var (elutionMin, elutionMax) = EstimatePeptideElutionWindow(scanIndex);
            double elutionSpan = elutionMax - elutionMin;
            Report($"  [Calibration] Elution window: {elutionMin:F2}–{elutionMax:F2} min (span={elutionSpan:F2} min)");

            // ── Sort all targets by iRT rank ───────────────────────────────────────
            var allTargetsSorted = new List<(int Idx, double LibRt)>();
            for (int i = 0; i < precursors.Count; i++)
            {
                var p = precursors[i];
                if (p.IsDecoy) continue;
                double libRt = p.IrtValue.HasValue ? p.IrtValue.Value
                             : p.RetentionTime.HasValue ? p.RetentionTime.Value
                             : double.NaN;
                if (!double.IsNaN(libRt))
                    allTargetsSorted.Add((i, libRt));
            }
            allTargetsSorted.Sort((a, b) => a.LibRt.CompareTo(b.LibRt));
            int totalTargets = allTargetsSorted.Count;

            Report($"  [Calibration] Sorted targets: {totalTargets:N0} with iRT values");

            if (totalTargets < MinAnchorCount)
            {
                Report("  [Calibration] Insufficient targets for center-outward bootstrap — falling back.");
                return LegacyProgressiveBootstrap(precursors, scanIndex, baseParameters, orchestrator, log);
            }
            var allBootstrapMedians = new List<(double Irt, double ApexRt)>();
            // ── Phase A: iterate center band (rank 40–60%) until slope converges ──
            // First pass uses rank-mapped flat windows (no model needed).
            // Subsequent passes use the current model to narrow windows around
            // predicted RT, improving slope accuracy before expanding outward.
            int centerLo = (int)(totalTargets * 0.40);
            int centerHi = (int)(totalTargets * 0.60);
            var centerIndices = new List<int>();
            for (int i = centerLo; i < centerHi; i++)
                centerIndices.Add(allTargetsSorted[i].Idx);

            double rtCenterLo = elutionMin + elutionSpan * 0.40;
            double rtCenterHi = elutionMin + elutionSpan * 0.60;
            double centerCellHalfWidth = Math.Max((rtCenterHi - rtCenterLo) / 5.0, MinWindowHalfWidthMinutes);

            var centerBand = BuildBandPrecursorList(precursors, centerIndices, includeDecoys: true);

            IRtCalibrationModel model0 = null;
            double prevSlopeA = double.NaN;
            List<DiaSearchResult> centerResults = null;

            for (int phaseAIter = 0; phaseAIter < 5; phaseAIter++)
            {
                DiaLibraryQueryGenerator.GenerationResult gen;
                double halfWidth;

                if (model0 == null)
                {
                    // First pass: rank-mapped flat windows, no model required
                    halfWidth = centerCellHalfWidth;
                    gen = GenerateWithRankMappedWindows(
                        centerBand, centerIndices, allTargetsSorted,
                        bandLo: centerLo, bandHi: centerHi,
                        rtBandLo: rtCenterLo, rtBandHi: rtCenterHi,
                        cellHalfWidth: halfWidth,
                        scanIndex, baseParameters);
                }
                else
                {
                    // Subsequent passes: model-predicted RT ± 2σ
                    halfWidth = Math.Max(model0.SigmaMinutes * 2.0, MinWindowHalfWidthMinutes);
                    gen = GenerateWithModelWindows(
                        centerBand, model0,
                        cellHalfWidth: halfWidth,
                        scanIndex, baseParameters);
                }

                extractSw.Restart();
                var extract = orchestrator.ExtractAll(gen.Queries);
                centerResults = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    centerBand, gen, extract, baseParameters, scanIndex);
                extractSw.Stop();

                fitSw.Restart();
                // Phase A iter 0 (first pass — rank-mapped, no model yet)
                var candidate = ComputeGridModel(centerResults, centerBand,
                    accumulator: allBootstrapMedians);
                fitSw.Stop();

                if (candidate == null)
                {
                    if (model0 == null)
                    {
                        Report($"  [Calibration] Phase A iter {phaseAIter} failed — falling back.");
                        return LegacyProgressiveBootstrap(precursors, scanIndex, baseParameters, orchestrator, log);
                    }
                    Report($"  [Calibration] Phase A iter {phaseAIter} grid failed — keeping previous model.");
                    break;
                }

                Report($"  [Calibration] Phase A iter {phaseAIter}: slope={candidate.Slope:F4}, intercept={candidate.Intercept:F3}, σ={candidate.SigmaMinutes:F3}, R²={candidate.RSquared:F4} (cell±{halfWidth:F2} min)");

                model0 = candidate;

                // Converged?
                if (!double.IsNaN(prevSlopeA) && Math.Abs(candidate.Slope - prevSlopeA) / prevSlopeA < 0.02)
                {
                    Report($"  [Calibration] Phase A converged at iter {phaseAIter}.");
                    break;
                }
                prevSlopeA = candidate.Slope;
            }

            Report($"  [Calibration] Phase A final: slope={model0.Slope:F4}, intercept={model0.Intercept:F3}, σ={model0.SigmaMinutes:F3}");


            // ── Phase B: mid bands (rank 20–40%, 60–80%) → model_1 ────────────────
            // Project Phase A model into B range — no blind grid search.
            // Use model-predicted windows and iterate until convergence.
            int midLoLo = (int)(totalTargets * 0.20);
            int midLoHi = centerLo;
            int midHiLo = centerHi;
            int midHiHi = (int)(totalTargets * 0.80);

            var midIndicesLo = new List<int>();
            var midIndicesHi = new List<int>();
            for (int i = midLoLo; i < midLoHi; i++) midIndicesLo.Add(allTargetsSorted[i].Idx);
            for (int i = midHiLo; i < midHiHi; i++) midIndicesHi.Add(allTargetsSorted[i].Idx);
            var midIndices = new List<int>(midIndicesLo.Count + midIndicesHi.Count);
            midIndices.AddRange(midIndicesLo);
            midIndices.AddRange(midIndicesHi);

            var midBand = BuildBandPrecursorList(precursors, midIndices, includeDecoys: true);

            IRtCalibrationModel model1 = model0;
            double prevSlopeB = double.NaN;
            List<DiaSearchResult> midResults = null;

            for (int phaseBIter = 0; phaseBIter < 5; phaseBIter++)
            {
                double halfWidth = Math.Max(model1.SigmaMinutes * 2.0, MinWindowHalfWidthMinutes);
                var gen = GenerateWithModelWindows(
                    midBand, model1,
                    cellHalfWidth: halfWidth,
                    scanIndex, baseParameters);

                extractSw.Restart();
                var extract = orchestrator.ExtractAll(gen.Queries);
                midResults = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    midBand, gen, extract, baseParameters, scanIndex);
                extractSw.Stop();

                var midBandLo = BuildBandPrecursorList(precursors, midIndicesLo, includeDecoys: false);
                var midBandHi = BuildBandPrecursorList(precursors, midIndicesHi, includeDecoys: false);
                var midLoKeys = new HashSet<(string, int)>(midBandLo.Select(p => (p.Sequence, p.ChargeState)));
                var midResultsLo = midResults.Where(r => midLoKeys.Contains((r.Sequence, r.ChargeState))).ToList();
                var midResultsHi = midResults.Where(r => !midLoKeys.Contains((r.Sequence, r.ChargeState))).ToList();

                fitSw.Restart();
                // Phase A iter 1+ (model-predicted windows)
                var candidate = ComputeGridModel(centerResults, centerBand,
                    accumulator: allBootstrapMedians);
                fitSw.Stop();

                if (candidate == null)
                {
                    Report($"  [Calibration] Phase B iter {phaseBIter} grid failed — keeping previous model.");
                    break;
                }

                Report($"  [Calibration] Phase B iter {phaseBIter}: slope={candidate.Slope:F4}, intercept={candidate.Intercept:F3}, σ={candidate.SigmaMinutes:F3}, R²={candidate.RSquared:F4} (cell±{halfWidth:F2} min)");

                model1 = candidate;

                if (!double.IsNaN(prevSlopeB) && Math.Abs(candidate.Slope - prevSlopeB) / prevSlopeB < 0.02)
                {
                    Report($"  [Calibration] Phase B converged at iter {phaseBIter}.");
                    break;
                }
                prevSlopeB = candidate.Slope;
            }

            Report($"  [Calibration] Phase B final: slope={model1.Slope:F4}, intercept={model1.Intercept:F3}, σ={model1.SigmaMinutes:F3}");

            // ── Phase C: outer bands (rank 0–20%, 80–100%) → currentModel ─────────
            // Project Phase B model into C range and iterate.
            int outerLoLo = 0;
            int outerLoHi = midLoLo;
            int outerHiLo = midHiHi;
            int outerHiHi = totalTargets;

            var outerIndicesLo = new List<int>();
            var outerIndicesHi = new List<int>();
            for (int i = outerLoLo; i < outerLoHi; i++) outerIndicesLo.Add(allTargetsSorted[i].Idx);
            for (int i = outerHiLo; i < outerHiHi; i++) outerIndicesHi.Add(allTargetsSorted[i].Idx);
            var outerIndices = new List<int>(outerIndicesLo.Count + outerIndicesHi.Count);
            outerIndices.AddRange(outerIndicesLo);
            outerIndices.AddRange(outerIndicesHi);

            var outerBand = BuildBandPrecursorList(precursors, outerIndices, includeDecoys: true);

            IRtCalibrationModel modelC = model1;
            double prevSlopeC = double.NaN;
            List<DiaSearchResult> outerResults = null;

            for (int phaseCIter = 0; phaseCIter < 5; phaseCIter++)
            {
                double halfWidth = Math.Max(modelC.SigmaMinutes * 2.0, MinWindowHalfWidthMinutes);
                var gen = GenerateWithModelWindows(
                    outerBand, modelC,
                    cellHalfWidth: halfWidth,
                    scanIndex, baseParameters);

                extractSw.Restart();
                var extract = orchestrator.ExtractAll(gen.Queries);
                outerResults = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    outerBand, gen, extract, baseParameters, scanIndex);
                extractSw.Stop();

                var outerBandLo = BuildBandPrecursorList(precursors, outerIndicesLo, includeDecoys: false);
                var outerBandHi = BuildBandPrecursorList(precursors, outerIndicesHi, includeDecoys: false);
                var outerLoKeys = new HashSet<(string, int)>(outerBandLo.Select(p => (p.Sequence, p.ChargeState)));
                var outerResultsLo = outerResults.Where(r => outerLoKeys.Contains((r.Sequence, r.ChargeState))).ToList();
                var outerResultsHi = outerResults.Where(r => !outerLoKeys.Contains((r.Sequence, r.ChargeState))).ToList();

                fitSw.Restart();
                var candidate = ComputeGridModelSplit(
    outerResultsLo, outerBandLo, outerResultsHi, outerBandHi,
    accumulator: allBootstrapMedians);
                fitSw.Stop();

                if (candidate == null)
                {
                    Report($"  [Calibration] Phase C iter {phaseCIter} grid failed — keeping previous model.");
                    break;
                }

                Report($"  [Calibration] Phase C iter {phaseCIter}: slope={candidate.Slope:F4}, intercept={candidate.Intercept:F3}, σ={candidate.SigmaMinutes:F3}, R²={candidate.RSquared:F4} (cell±{halfWidth:F2} min)");

                modelC = candidate;

                if (!double.IsNaN(prevSlopeC) && Math.Abs(candidate.Slope - prevSlopeC) / prevSlopeC < 0.02)
                {
                    Report($"  [Calibration] Phase C converged at iter {phaseCIter}.");
                    break;
                }
                prevSlopeC = candidate.Slope;
            }

            // ── Select final bootstrap model and cap sigma ────────────────────────
            const double BootstrapSigmaCap = 0.3;

            IRtCalibrationModel bootstrapFinal = modelC;

            if (EnableNonLinearModelSelection && allBootstrapMedians.Count >= 50)
            {
                var libRts = allBootstrapMedians.Select(m => m.Irt).ToArray();
                var apexRts = allBootstrapMedians.Select(m => m.ApexRt).ToArray();

                var pwModel = PiecewiseLinearRtModel.Fit(libRts, apexRts, numSegments: 8);

                if (pwModel != null && pwModel.SigmaMinutes < modelC.SigmaMinutes
                                    && pwModel.RSquared >= 0.995)
                {
                    Report($"  [Calibration] Piecewise model selected: σ={pwModel.SigmaMinutes:F3}, R²={pwModel.RSquared:F4}, segments={pwModel.NumSegments}");
                    bootstrapFinal = pwModel;
                }
                else
                {
                    Report($"  [Calibration] Linear model retained: σ={modelC.SigmaMinutes:F3}, R²={modelC.RSquared:F4}");
                }
            }

            if (bootstrapFinal.SigmaMinutes > BootstrapSigmaCap)
            {
                var capped = new RtCalibrationModel(
                    slope: bootstrapFinal.Slope,
                    intercept: bootstrapFinal.Intercept,
                    sigmaMinutes: BootstrapSigmaCap,
                    rSquared: bootstrapFinal.RSquared,
                    anchorCount: bootstrapFinal.AnchorCount);
                currentModel = new LinearRtModelWrapper(capped);
                Report($"  [Calibration] Bootstrap complete: slope={bootstrapFinal.Slope:F4}, intercept={bootstrapFinal.Intercept:F3}, σ={bootstrapFinal.SigmaMinutes:F3}→{BootstrapSigmaCap:F3} (capped), R²={bootstrapFinal.RSquared:F4}");
            }
            else
            {
                currentModel = bootstrapFinal;
                Report($"  [Calibration] Bootstrap complete: slope={currentModel.Slope:F4}, intercept={currentModel.Intercept:F3}, σ={currentModel.SigmaMinutes:F3}, R²={currentModel.RSquared:F4}");
            }


            // Combine all results for refinement iterations
            currentResults = new List<DiaSearchResult>(centerResults.Count + midResults.Count + outerResults.Count);
            currentResults.AddRange(centerResults);
            currentResults.AddRange(midResults);
            currentResults.AddRange(outerResults);

            previousSigma = currentModel.SigmaMinutes;
            previousSlope = currentModel.Slope;

            double hw0 = Math.Max(currentModel.SigmaMinutes * SigmaMultiplier, MinWindowHalfWidthMinutes);
            log.Add(CreateLogEntry(0, currentResults.Count, currentModel, hw0,
                extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, "CenterOutward"));

            double bootstrapWindowCap = Math.Max(currentModel.SigmaMinutes * BootstrapSigmaMultiplier, MinRefinementWindowMinutes);
            return RunRefinementIterations(
                precursors, scanIndex, baseParameters, orchestrator,
                currentModel, currentResults, previousSigma, previousSlope,
                log, startIteration: 1,
                maxWindowHalfWidth: bootstrapWindowCap);
        }

        // ── Offset Detection Helpers ───────────────────────────────────────

        private static double FindOffsetMode(List<double> offsets, double bandwidth)
        {
            if (offsets.Count == 0) return double.NaN;

            // Clamp grid to IQR ± 3×bandwidth to avoid outliers diluting resolution.
            // A handful of extreme false matches can drag min/max far from the true mode,
            // spreading 500 grid points over a range where we have no signal.
            var sorted = offsets.OrderBy(o => o).ToList();
            int n = sorted.Count;
            double q1 = sorted[n / 4];
            double q3 = sorted[3 * n / 4];
            double iqr = q3 - q1;
            double gridMin = q1 - Math.Max(iqr * 1.5, 3.0 * bandwidth);
            double gridMax = q3 + Math.Max(iqr * 1.5, 3.0 * bandwidth);
            if (gridMax - gridMin < 1e-6) return sorted[n / 2];

            const int gridPoints = 500;
            double step = (gridMax - gridMin) / (gridPoints - 1);
            double bestDensity = double.MinValue, bestX = double.NaN;
            for (int g = 0; g < gridPoints; g++)
            {
                double x = gridMin + g * step;
                double density = 0.0;
                foreach (double o in offsets)
                {
                    double d = (x - o) / bandwidth;
                    density += Math.Exp(-0.5 * d * d);
                }
                if (density > bestDensity) { bestDensity = density; bestX = x; }
            }
            return bestX;
        }

        private static List<LibraryPrecursorInput> BuildSubsample(
            IList<LibraryPrecursorInput> precursors, int maxSize)
        {
            // Collect targets only — decoys have no reliable RT for offset detection.
            // Sort by RetentionTime/IrtValue so reservoir sampling draws uniformly
            // across the RT range, not from whatever order the library was written.
            var targets = new List<LibraryPrecursorInput>(precursors.Count / 2 + 1);
            foreach (var p in precursors)
                if (!p.IsDecoy) targets.Add(p);

            targets.Sort((a, b) =>
            {
                double rtA = a.IrtValue ?? a.RetentionTime ?? 0;
                double rtB = b.IrtValue ?? b.RetentionTime ?? 0;
                return rtA.CompareTo(rtB);
            });

            if (targets.Count <= maxSize) return targets;

            // Evenly-spaced stride gives better RT coverage than random reservoir
            // when the list is already RT-sorted.
            double stride = (double)targets.Count / maxSize;
            var result = new List<LibraryPrecursorInput>(maxSize);
            for (int i = 0; i < maxSize; i++)
                result.Add(targets[(int)(i * stride)]);
            return result;
        }

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

            // Skip decoys: offset bootstrap is used only to fit a calibration model,
            // not to score decoys. Including decoys wastes extraction time and inflates
            // the result list with hits that can't contribute a reliable LibraryRt offset.
            int totalQueryCount = 0, skippedNoWindow = 0, skippedNoFragments = 0;
            for (int i = 0; i < precursors.Count; i++)
            {
                var p = precursors[i];
                if (p.IsDecoy) continue;
                int windowId = scanIndex.FindWindowForPrecursorMz(p.PrecursorMz);
                if (windowId < 0) { skippedNoWindow++; continue; }
                if (p.FragmentCount == 0) { skippedNoFragments++; continue; }
                totalQueryCount += p.FragmentCount;
            }

            var queries = new FragmentQuery[totalQueryCount];
            var groups = new List<DiaLibraryQueryGenerator.PrecursorQueryGroup>(
                precursors.Count / 2);

            int queryIndex = 0;
            for (int i = 0; i < precursors.Count; i++)
            {
                var p = precursors[i];
                if (p.IsDecoy) continue;
                int windowId = scanIndex.FindWindowForPrecursorMz(p.PrecursorMz);
                if (windowId < 0 || p.FragmentCount == 0) continue;

                double libRt = p.IrtValue.HasValue ? p.IrtValue.Value
                             : p.RetentionTime.HasValue ? p.RetentionTime.Value
                             : double.NaN;

                float rtMin, rtMax;
                if (!double.IsNaN(libRt))
                {
                    float center = (float)(libRt + offset);
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
                int startIteration,
                double maxWindowHalfWidth = double.MaxValue)
        {
            for (int iteration = startIteration; iteration < MaxIterations; iteration++)
            {
                Console.WriteLine($"  [Calibration] Iteration {iteration} starting...");
                var extractSw = new Stopwatch();
                var fitSw = new Stopwatch();
                var selectSw = new Stopwatch();

                extractSw.Start();

                // Window may only shrink from the bootstrap window — never grow beyond it.
                // Using SigmaMultiplier × σ without a cap caused the window to expand on
                // the first refinement iteration (bootstrap used 3σ, refinement uses 4σ),
                // pulling in more false matches, inflating σ further — a divergence loop.
                double rawWindow = Math.Max(
                    currentModel.SigmaMinutes * SigmaMultiplier,
                    MinRefinementWindowMinutes);
                double windowHalfWidth = Math.Min(rawWindow, maxWindowHalfWidth);

                var genResult = GenerateCalibratedWithExplicitWindow(
                    precursors, scanIndex, baseParameters, currentModel, windowHalfWidth);
                var extractionResult = orchestrator.ExtractAll(genResult.Queries);
                var results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    precursors, genResult, extractionResult, baseParameters, scanIndex);
                extractSw.Stop();

                Console.WriteLine($"  [Calibration] Iteration {iteration}: extraction={extractSw.ElapsedMilliseconds}ms, results={results.Count:N0}");
                ReportTdRatio($"Iter {iteration}", results, RefinedApexScoreThreshold);

                selectSw.Start();
                var anchors = SelectAnchors(results, precursors, genResult.PrecursorGroups,
                    RefinedApexScoreThreshold, int.MaxValue, requireMinFragDetRate: true,
                    logger: Report);
                selectSw.Stop();

                Console.WriteLine($"  [Calibration] Iteration {iteration}: {anchors.Count} anchors (threshold={RefinedApexScoreThreshold:F2})");

                if (anchors.Count < MinAnchorCount)
                {
                    log.Add(CreateLogEntry(iteration, anchors.Count, currentModel,
                        Math.Min(Math.Max(currentModel.SigmaMinutes * SigmaMultiplier, MinRefinementWindowMinutes), maxWindowHalfWidth),
                        extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, "InsufficientAnchors"));
                    break;
                }

                fitSw.Start();
                var libraryRts = anchors.Select(a => a.LibraryRt).ToArray();
                var observedRts = anchors.Select(a => a.ObservedRt).ToArray();

                // RANSAC threshold: start tight (0.5 min covers ~2-3× true σ of 0.15–0.25 min).
                // The old 2.0×σ threshold was using the inflated σ from the wide inlier set,
                // keeping false matches in the fit indefinitely.
                double inlierThreshold = Math.Max(currentModel.SigmaMinutes * 2.0, 0.3);
                var linearModel = FitLinearWithRansac(libraryRts, observedRts, inlierThreshold);

                if (linearModel == null || !linearModel.IsReliable)
                {
                    log.Add(CreateLogEntry(iteration, anchors.Count, currentModel,
                        Math.Min(Math.Max(currentModel.SigmaMinutes * SigmaMultiplier, MinRefinementWindowMinutes), maxWindowHalfWidth),
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
                        Math.Min(Math.Max(newSigma * SigmaMultiplier, MinRefinementWindowMinutes), maxWindowHalfWidth),
                        extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, modelTypeLabel + "_Diverged_Reverted"));
                    break;
                }

                // Slope stability: each iteration should only refine, not re-estimate wildly.
                // If slope moves >10% relative to previous, the RANSAC found a bad consensus.
                if (!double.IsNaN(previousSlope) && Math.Abs(previousSlope) > 1e-6)
                {
                    double relSlopeChange = Math.Abs(newSlope - previousSlope) / Math.Abs(previousSlope);
                    if (relSlopeChange > 0.10)
                    {
                        Report($"  [Calibration] Iteration {iteration}: slope changed {relSlopeChange:P0} ({previousSlope:F4}→{newSlope:F4}), reverting.");
                        log.Add(CreateLogEntry(iteration, anchors.Count, currentModel,
                            Math.Min(Math.Max(currentModel.SigmaMinutes * SigmaMultiplier, MinRefinementWindowMinutes), maxWindowHalfWidth),
                            extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, modelTypeLabel + "_SlopeDiverged_Reverted"));
                        break;
                    }
                }

                currentModel = bestModel;
                currentResults = results;

                double hw = Math.Min(
                    Math.Min(Math.Max(newSigma * SigmaMultiplier, MinRefinementWindowMinutes), maxWindowHalfWidth),
                    maxWindowHalfWidth);
                log.Add(CreateLogEntry(iteration, anchors.Count, bestModel, hw,
                    extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, modelTypeLabel));

                if (iteration > startIteration && previousSigma > 0 && !double.IsInfinity(previousSigma))
                {
                    double relativeDelta = Math.Abs(newSigma - previousSigma) / previousSigma;
                    if (relativeDelta < ConvergenceThreshold)
                    {
                        Report($"  [Calibration] Converged at iteration {iteration} (Δσ/σ = {relativeDelta:F4} < {ConvergenceThreshold})");
                        break;
                    }
                }

                previousSigma = newSigma;
                previousSlope = newSlope;
                // Ratchet the window cap down as σ improves — never let it grow back up.
                // This ensures refinement always narrows toward the true RT distribution.
                double nextRawWindow = Math.Max(newSigma * SigmaMultiplier, MinRefinementWindowMinutes);
                if (nextRawWindow < maxWindowHalfWidth)
                    maxWindowHalfWidth = nextRawWindow;
            }

            if (currentModel == null)
                return (null, null, log, currentResults);

            return (currentModel.ToRtCalibrationModel(), currentModel, log, currentResults);
        }

        // ── Phase 23 Bootstrap Helpers ─────────────────────────────────────

        /// <summary>
        /// Estimates the peptide elution window from the MS2 TIC profile.
        /// Sums fragment intensities per scan, applies a 5-scan median smooth,
        /// and returns the RT range where the smoothed TIC exceeds TicThresholdFraction
        /// of the peak TIC. Falls back to the global RT range if detection fails.
        /// </summary>
        /// <summary>
        /// Estimates the peptide elution window using cumulative TIC.
        /// Returns the RT where 10% of cumulative TIC has passed (gradient start)
        /// and the RT where 90% of cumulative TIC has passed (gradient end).
        /// This brackets the actual peptide elution regardless of gradient shape.
        /// </summary>
        private (double Min, double Max) EstimatePeptideElutionWindow(DiaScanIndex scanIndex)
        {
            // Use the actual file RT range trimmed by a fixed margin.
            //
            // Why not TIC percentiles: the MS2 TIC is dominated by the dense gradient
            // midpoint. Even at 2/98 percentile, the window ends up being ~60% of the
            // true gradient span, which compresses all rank-mapped extraction windows
            // and causes inverted T/D ratios throughout the bootstrap.
            //
            // The file RT range already excludes pre-column equilibration dead time
            // (vendors set scan start/end to the gradient boundaries). We apply a
            // 2-minute fixed trim on each side as a conservative guard against the
            // very first/last few scans being noisy survey spectra.

            double rtMin = scanIndex.GetGlobalRtMin();
            double rtMax = scanIndex.GetGlobalRtMax();
            double span = rtMax - rtMin;

            const double TrimMinutes = 2.0;

            // Only trim if the run is long enough that trimming makes sense
            if (span > TrimMinutes * 4)
            {
                rtMin += TrimMinutes;
                rtMax -= TrimMinutes;
            }

            return (rtMin, rtMax);
        }
        /// <summary>
        /// Builds a precursor list for a single iRT band.
        /// Always includes all decoys so FDR scoring works correctly.
        /// targetIndices are original indices into the full precursors list.
        /// </summary>
        private static List<LibraryPrecursorInput> BuildBandPrecursorList(
            IList<LibraryPrecursorInput> all,
            IList<int> targetIndices,
            bool includeDecoys)
        {
            var result = new List<LibraryPrecursorInput>(
                targetIndices.Count * (includeDecoys ? 2 : 1));

            foreach (int idx in targetIndices)
                result.Add(all[idx]);

            if (!includeDecoys)
                return result;

            // Find where decoys start in the combined list (targets come first)
            int firstDecoyIdx = -1;
            for (int i = 0; i < all.Count; i++)
            {
                if (all[i].IsDecoy) { firstDecoyIdx = i; break; }
            }

            if (firstDecoyIdx < 0)
                return result; // no decoys in library

            // For each target in this band, include its paired decoy.
            // CRITICAL: decoys must inherit the target's IrtValue/RetentionTime, NOT
            // their own — decoy sequences are scrambled so their iRT predictions are
            // unrelated to the target's observed RT. Using the decoy's own RT would
            // place it in a completely different RT space, making grid T/D comparisons
            // meaningless and causing RANSAC to fail.
            foreach (int targetIdx in targetIndices)
            {
                var target = all[targetIdx];

                // Paired decoy is at offset targetIdx within the decoy block
                int decoyIdx = firstDecoyIdx + targetIdx;
                if (decoyIdx >= all.Count || !all[decoyIdx].IsDecoy)
                    continue;

                var decoy = all[decoyIdx];

                // Rebuild the decoy with the target's RT metadata so it gets
                // searched in the correct RT window
                result.Add(new LibraryPrecursorInput(
                    sequence: decoy.Sequence,
                    precursorMz: decoy.PrecursorMz,
                    chargeState: decoy.ChargeState,
                    retentionTime: target.RetentionTime,   // ← target's RT
                    isDecoy: true,
                    fragmentMzs: decoy.FragmentMzs,
                    fragmentIntensities: decoy.FragmentIntensities,
                    irtValue: target.IrtValue));            // ← target's iRT
            }

            return result;
        }

        /// <summary>
        /// Generates queries with a fixed [rtMin, rtMax] for every precursor.
        /// Used for Phase A where no calibration model exists yet.
        /// Every precursor in the list receives the same RT window.
        /// </summary>
        private static DiaLibraryQueryGenerator.GenerationResult GenerateWithFlatWindowForAll(
            IList<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            DiaSearchParameters parameters,
            float rtMin,
            float rtMax)
        {
            float ppmTolerance = parameters.PpmTolerance;
            float globalRtMin = scanIndex.GetGlobalRtMin();
            float globalRtMax = scanIndex.GetGlobalRtMax();

            // Clamp to actual run range
            float clampedMin = Math.Max(rtMin, globalRtMin);
            float clampedMax = Math.Min(rtMax, globalRtMax);
            if (clampedMin >= clampedMax) { clampedMin = globalRtMin; clampedMax = globalRtMax; }

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

                int groupOffset = queryIndex;
                for (int f = 0; f < p.FragmentCount; f++)
                {
                    queries[queryIndex] = new FragmentQuery(
                        targetMz: p.FragmentMzs[f],
                        tolerancePpm: ppmTolerance,
                        rtMin: clampedMin,
                        rtMax: clampedMax,
                        windowId: windowId,
                        queryId: queryIndex);
                    queryIndex++;
                }

                groups.Add(new DiaLibraryQueryGenerator.PrecursorQueryGroup(
                    inputIndex: i,
                    queryOffset: groupOffset,
                    queryCount: p.FragmentCount,
                    windowId: windowId,
                    rtMin: clampedMin,
                    rtMax: clampedMax));
            }

            return new DiaLibraryQueryGenerator.GenerationResult(
                queries, groups.ToArray(), skippedNoWindow, skippedNoFragments);
        }

        /// <summary>
        /// Generates queries for the center-band targets where each precursor's RT window
        /// is determined by its iRT rank position within the band, mapped linearly onto
        /// the corresponding RT sub-band.
        ///
        /// For a precursor at rank fraction f within the band (0=lowest iRT, 1=highest),
        /// its expected RT is: rtBandLo + f * rtBandSpan.
        /// Its window is that predicted RT ± cellHalfWidth.
        ///
        /// This gives every peptide a narrow, correctly-positioned RT window without
        /// requiring a calibration model. It works because the mapping from iRT rank
        /// to RT rank is monotone regardless of the unknown slope.
        ///
        /// targetIndices: original indices into 'all' for the band targets, in iRT-sorted order.
        /// bandLo/bandHi: the rank range within the global sorted targets list.
        /// </summary>
        private static DiaLibraryQueryGenerator.GenerationResult GenerateWithRankMappedWindows(
            IList<LibraryPrecursorInput> bandPrecursors,
            IList<int> targetIndices,
            IList<(int Idx, double LibRt)> allTargetsSorted,
            int bandLo,
            int bandHi,
            double rtBandLo,
            double rtBandHi,
            double cellHalfWidth,
            DiaScanIndex scanIndex,
            DiaSearchParameters parameters)
        {
            float ppmTolerance = parameters.PpmTolerance;
            float globalRtMin = scanIndex.GetGlobalRtMin();
            float globalRtMax = scanIndex.GetGlobalRtMax();

            double rtBandSpan = rtBandHi - rtBandLo;
            int bandCount = bandHi - bandLo;

            // Build a lookup: original precursor index → rank-fraction within band
            // targetIndices[j] == allTargetsSorted[bandLo + j].Idx
            var rankFraction = new Dictionary<int, double>(targetIndices.Count);
            for (int j = 0; j < targetIndices.Count; j++)
            {
                double f = bandCount > 1 ? (double)j / (bandCount - 1) : 0.5;
                rankFraction[targetIndices[j]] = f;
            }

            // Count queries
            int totalQueryCount = 0, skippedNoWindow = 0, skippedNoFragments = 0;
            for (int i = 0; i < bandPrecursors.Count; i++)
            {
                var p = bandPrecursors[i];
                int windowId = scanIndex.FindWindowForPrecursorMz(p.PrecursorMz);
                if (windowId < 0) { skippedNoWindow++; continue; }
                if (p.FragmentCount == 0) { skippedNoFragments++; continue; }
                totalQueryCount += p.FragmentCount;
            }

            var queries = new FragmentQuery[totalQueryCount];
            var groups = new List<DiaLibraryQueryGenerator.PrecursorQueryGroup>(
                bandPrecursors.Count - skippedNoWindow - skippedNoFragments);

            // We need to map bandPrecursors[i] back to its original index.
            // BuildBandPrecursorList puts targets first (in targetIndices order), then decoys.
            // So bandPrecursors[j] for j < targetIndices.Count corresponds to targetIndices[j].
            int queryIndex = 0;
            for (int i = 0; i < bandPrecursors.Count; i++)
            {
                var p = bandPrecursors[i];
                int windowId = scanIndex.FindWindowForPrecursorMz(p.PrecursorMz);
                if (windowId < 0 || p.FragmentCount == 0)
                    continue;

                float rtMin, rtMax;

                // Is this a target in the band with a known rank fraction?
                // Targets come first in the list (indices 0..targetIndices.Count-1).
                // Decoys get the SAME narrow window as their corresponding rank-mapped target
                // (matched by position in the list after the target block). This prevents
                // the full-band decoy flood that destroys the T/D signal in the grid cells.
                double rankF;
                if (i < targetIndices.Count && rankFraction.TryGetValue(targetIndices[i], out rankF))
                {
                    // This is a band target — use its rank-mapped window
                    double predictedRt = rtBandLo + rankF * rtBandSpan;
                    rtMin = (float)Math.Max(predictedRt - cellHalfWidth, globalRtMin);
                    rtMax = (float)Math.Min(predictedRt + cellHalfWidth, globalRtMax);
                }
                else
                {
                    // Decoy: pair with the matching band target by position offset.
                    // i - targetIndices.Count gives the decoy's position in the appended decoy block.
                    // We cycle through rank fractions so each decoy is searched in a specific cell,
                    // not the whole band. This keeps decoy counts proportional to cell size (fair comparison).
                    int decoyPos = i - targetIndices.Count;
                    double decoyF = bandCount > 1 ? (double)(decoyPos % bandCount) / (bandCount - 1) : 0.5;
                    double decoyPredictedRt = rtBandLo + decoyF * rtBandSpan;
                    rtMin = (float)Math.Max(decoyPredictedRt - cellHalfWidth, globalRtMin);
                    rtMax = (float)Math.Min(decoyPredictedRt + cellHalfWidth, globalRtMax);
                }

                int groupOffset = queryIndex;
                for (int f2 = 0; f2 < p.FragmentCount; f2++)
                {
                    queries[queryIndex] = new FragmentQuery(
                        targetMz: p.FragmentMzs[f2],
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
        private static DiaLibraryQueryGenerator.GenerationResult GenerateWithModelWindows(
    IList<LibraryPrecursorInput> bandPrecursors,
    IRtCalibrationModel model,
    double cellHalfWidth,
    DiaScanIndex scanIndex,
    DiaSearchParameters parameters)
        {
            float ppmTolerance = parameters.PpmTolerance;
            float globalRtMin = scanIndex.GetGlobalRtMin();
            float globalRtMax = scanIndex.GetGlobalRtMax();

            int totalQueryCount = 0, skippedNoWindow = 0, skippedNoFragments = 0;
            for (int i = 0; i < bandPrecursors.Count; i++)
            {
                var p = bandPrecursors[i];
                int windowId = scanIndex.FindWindowForPrecursorMz(p.PrecursorMz);
                if (windowId < 0) { skippedNoWindow++; continue; }
                if (p.FragmentCount == 0) { skippedNoFragments++; continue; }
                totalQueryCount += p.FragmentCount;
            }

            var queries = new FragmentQuery[totalQueryCount];
            var groups = new List<DiaLibraryQueryGenerator.PrecursorQueryGroup>(
                bandPrecursors.Count - skippedNoWindow - skippedNoFragments);

            int queryIndex = 0;
            for (int i = 0; i < bandPrecursors.Count; i++)
            {
                var p = bandPrecursors[i];
                int windowId = scanIndex.FindWindowForPrecursorMz(p.PrecursorMz);
                if (windowId < 0 || p.FragmentCount == 0)
                    continue;

                // Get library RT for this precursor
                double libRt = p.IrtValue.HasValue ? p.IrtValue.Value
                             : p.RetentionTime.HasValue ? p.RetentionTime.Value
                             : double.NaN;

                float rtMin, rtMax;
                if (!double.IsNaN(libRt))
                {
                    double predictedRt = model.ToMinutes(libRt);
                    rtMin = (float)Math.Max(predictedRt - cellHalfWidth, globalRtMin);
                    rtMax = (float)Math.Min(predictedRt + cellHalfWidth, globalRtMax);
                }
                else
                {
                    // No iRT — search full elution range as fallback
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
        /// Fits an initial calibration model by scoring a T/D grid.
        ///
        /// Divides [rtBandLo, rtBandHi] and [irtBandLo, irtBandHi] each into numSegments
        /// cells, forming an (numSegments × numSegments) grid. For each cell, counts
        /// targets and decoys whose observed apex RT and library RT fall in that cell.
        ///
        /// Cells with T/D ratio well above 1:1 are "signal" cells. These form a diagonal
        /// whose slope defines the iRT→RT relationship for this band. A weighted OLS fit
        /// through the signal-cell centers gives the initial model.
        ///
        /// Returns null if fewer than 2 signal cells are found.
        /// </summary>
        private IRtCalibrationModel ComputeGridModel(
            List<DiaSearchResult> results,
            IList<LibraryPrecursorInput> bandPrecursors,
            int targetsPerBin = 50,
            List<(double Irt, double ApexRt)> accumulator = null)
        {
            // Build lookup: (Sequence, ChargeState) → library iRT from precursor's own fields
            var irtByKey = new Dictionary<(string, int), double>(bandPrecursors.Count);
            foreach (var p in bandPrecursors)
            {
                double libRt = p.IrtValue.HasValue ? p.IrtValue.Value
                             : p.RetentionTime.HasValue ? p.RetentionTime.Value
                             : double.NaN;
                if (!double.IsNaN(libRt))
                    irtByKey[(p.Sequence, p.ChargeState)] = libRt;
            }

            // Build lookup: (Sequence, ChargeState) → best apex RT among target results
            var apexByKey = new Dictionary<(string, int), (double ApexRt, float ApexScore)>();
            foreach (var r in results)
            {
                if (r.IsDecoy) continue;
                if (float.IsNaN(r.ObservedApexRt)) continue;
                var key = (r.Sequence, r.ChargeState);
                if (!apexByKey.TryGetValue(key, out var existing)
                    || r.ApexScore > existing.ApexScore)
                    apexByKey[key] = (r.ObservedApexRt, r.ApexScore);
            }

            // Only include non-decoy precursors that have iRT values.
            // bandPrecursors is already in iRT-sorted order.
            var targetPrecursors = new List<LibraryPrecursorInput>();
            foreach (var p in bandPrecursors)
                if (!p.IsDecoy && irtByKey.ContainsKey((p.Sequence, p.ChargeState)))
                    targetPrecursors.Add(p);

            int numBins = Math.Max(3, targetPrecursors.Count / targetsPerBin);
            int binSize = targetPrecursors.Count / numBins;

            var binMedianIrt = new List<double>();
            var binMedianRt = new List<double>();

            for (int b = 0; b < numBins; b++)
            {
                int lo = b * binSize;
                int hi = (b == numBins - 1) ? targetPrecursors.Count : lo + binSize;

                var binIrts = new List<double>();
                var binApexRts = new List<double>();

                for (int i = lo; i < hi; i++)
                {
                    var p = targetPrecursors[i];
                    var aKey = (p.Sequence, p.ChargeState);
                    if (!irtByKey.TryGetValue(aKey, out double irt)) continue;
                    if (!apexByKey.TryGetValue(aKey, out var apex)) continue;
                    binIrts.Add(irt);
                    binApexRts.Add(apex.ApexRt);
                }

                if (binIrts.Count < 5) continue;

                // Take median of iRT and median of apex RT independently
                binIrts.Sort();
                binApexRts.Sort();
                int mid = binIrts.Count / 2;
                binMedianIrt.Add(binIrts[mid]);
                binMedianRt.Add(binApexRts[mid]);
            }
            if (accumulator != null)
            {
                for (int i = 0; i < binMedianIrt.Count; i++)
                    accumulator.Add((binMedianIrt[i], binMedianRt[i]));
            }
            Report($"  [Calibration] Apex grid: {numBins} bins, {binMedianIrt.Count} valid (targetsPerBin={targetsPerBin})");

            // Diagnostic: log first 5 and last 5 bin medians to verify grid anchoring
            int logCount = Math.Min(5, binMedianIrt.Count);
            var diagLines = new System.Text.StringBuilder();
            diagLines.Append("  [Calibration] Grid bin medians (first/last 5): ");
            for (int i = 0; i < logCount; i++)
                diagLines.Append($"[{binMedianIrt[i]:F1}→{binMedianRt[i]:F2}] ");
            if (binMedianIrt.Count > 10)
                diagLines.Append("... ");
            for (int i = Math.Max(logCount, binMedianIrt.Count - 5); i < binMedianIrt.Count; i++)
                diagLines.Append($"[{binMedianIrt[i]:F1}→{binMedianRt[i]:F2}] ");
            Report(diagLines.ToString());

            if (binMedianIrt.Count < 3)
                return null;

            // OLS through bin medians
            double n = binMedianIrt.Count;
            double sumX = 0, sumY = 0, sumXX = 0, sumXY = 0;
            for (int i = 0; i < binMedianIrt.Count; i++)
            {
                sumX += binMedianIrt[i];
                sumY += binMedianRt[i];
                sumXX += binMedianIrt[i] * binMedianIrt[i];
                sumXY += binMedianIrt[i] * binMedianRt[i];
            }

            double denom = n * sumXX - sumX * sumX;
            if (Math.Abs(denom) < 1e-12) return null;

            double slope = (n * sumXY - sumX * sumY) / denom;
            double intercept = (sumY - slope * sumX) / n;

            if (slope <= 0) return null;

            double ssRes = 0, meanRt = sumY / n, ssTot = 0;
            for (int i = 0; i < binMedianIrt.Count; i++)
            {
                double pred = slope * binMedianIrt[i] + intercept;
                double res = binMedianRt[i] - pred;
                ssRes += res * res;
                ssTot += (binMedianRt[i] - meanRt) * (binMedianRt[i] - meanRt);
            }
            double sigma = binMedianIrt.Count > 2
                ? Math.Sqrt(ssRes / (binMedianIrt.Count - 2))
                : 0.5;
            double rSquared = ssTot > 0 ? Math.Max(0, 1.0 - ssRes / ssTot) : 0;

            Report($"  [Calibration] Apex grid model: slope={slope:F4}, intercept={intercept:F3}, σ={sigma:F3}, R²={rSquared:F4}");

            var syntheticModel = new RtCalibrationModel(
                slope: slope,
                intercept: intercept,
                sigmaMinutes: Math.Max(sigma, MinWindowHalfWidthMinutes),
                rSquared: rSquared,
                anchorCount: binMedianIrt.Count);

            return new LinearRtModelWrapper(syntheticModel);
        }
        private IRtCalibrationModel ComputeGridModelSplit(
    List<DiaSearchResult> resultsLo, IList<LibraryPrecursorInput> bandLo,
    List<DiaSearchResult> resultsHi, IList<LibraryPrecursorInput> bandHi,
    int targetsPerBin = 50,
    List<(double Irt, double ApexRt)> accumulator = null)
        {
            var allMedians = new List<(double Irt, double ApexRt)>();
            CollectBinMedians(resultsLo, bandLo, targetsPerBin, allMedians);
            CollectBinMedians(resultsHi, bandHi, targetsPerBin, allMedians);
            accumulator?.AddRange(allMedians);

            Report($"  [Calibration] Apex grid (split): {allMedians.Count} bin medians total");
            if (allMedians.Count < 3) return null;

            allMedians.Sort((a, b) => a.Irt.CompareTo(b.Irt));
            double n = allMedians.Count;
            double sumX = 0, sumY = 0, sumXX = 0, sumXY = 0;
            foreach (var (irt, rt) in allMedians)
            {
                sumX += irt; sumY += rt;
                sumXX += irt * irt; sumXY += irt * rt;
            }
            double denom = n * sumXX - sumX * sumX;
            if (Math.Abs(denom) < 1e-12) return null;
            double slope = (n * sumXY - sumX * sumY) / denom;
            double intercept = (sumY - slope * sumX) / n;
            if (slope <= 0) return null;

            double ssRes = 0, meanRt = sumY / n, ssTot = 0;
            foreach (var (irt, rt) in allMedians)
            {
                double pred = slope * irt + intercept;
                double res = rt - pred;
                ssRes += res * res;
                ssTot += (rt - meanRt) * (rt - meanRt);
            }
            double sigma = allMedians.Count > 2
                ? Math.Sqrt(ssRes / (allMedians.Count - 2)) : 0.5;
            double rSquared = ssTot > 0 ? Math.Max(0, 1.0 - ssRes / ssTot) : 0;

            Report($"  [Calibration] Apex grid model: slope={slope:F4}, intercept={intercept:F3}, σ={sigma:F3}, R²={rSquared:F4}");
            return new LinearRtModelWrapper(new RtCalibrationModel(
                slope: slope, intercept: intercept,
                sigmaMinutes: Math.Max(sigma, MinWindowHalfWidthMinutes),
                rSquared: rSquared, anchorCount: allMedians.Count));
        }

        private void CollectBinMedians(
    List<DiaSearchResult> results,
    IList<LibraryPrecursorInput> bandPrecursors,
    int targetsPerBin,
    List<(double Irt, double ApexRt)> output)
        {
            var irtByKey = new Dictionary<(string, int), double>(bandPrecursors.Count);
            foreach (var p in bandPrecursors)
            {
                double libRt = p.IrtValue.HasValue ? p.IrtValue.Value
                             : p.RetentionTime.HasValue ? p.RetentionTime.Value
                             : double.NaN;
                if (!double.IsNaN(libRt))
                    irtByKey[(p.Sequence, p.ChargeState)] = libRt;
            }
            var apexByKey = new Dictionary<(string, int), (double ApexRt, float ApexScore)>();
            foreach (var r in results)
            {
                if (r.IsDecoy || float.IsNaN(r.ObservedApexRt)) continue;
                var key = (r.Sequence, r.ChargeState);
                if (!apexByKey.TryGetValue(key, out var existing) || r.ApexScore > existing.ApexScore)
                    apexByKey[key] = (r.ObservedApexRt, r.ApexScore);
            }
            var targetPrecursors = new List<LibraryPrecursorInput>();
            foreach (var p in bandPrecursors)
                if (!p.IsDecoy && irtByKey.ContainsKey((p.Sequence, p.ChargeState)))
                    targetPrecursors.Add(p);

            int numBins = Math.Max(3, targetPrecursors.Count / targetsPerBin);
            int binSize = targetPrecursors.Count / numBins;

            for (int b = 0; b < numBins; b++)
            {
                int lo = b * binSize;
                int hi = (b == numBins - 1) ? targetPrecursors.Count : lo + binSize;

                var binIrts = new List<double>();
                var binApexRts = new List<double>();

                for (int i = lo; i < hi; i++)
                {
                    var p = targetPrecursors[i];
                    var aKey = (p.Sequence, p.ChargeState);
                    if (!irtByKey.TryGetValue(aKey, out double irt)) continue;
                    if (!apexByKey.TryGetValue(aKey, out var apex)) continue;
                    binIrts.Add(irt);
                    binApexRts.Add(apex.ApexRt);
                }

                if (binIrts.Count < 5) continue;

                // Independent medians — robust to outliers on either axis
                binIrts.Sort();
                binApexRts.Sort();
                int mid = binIrts.Count / 2;
                output.Add((binIrts[mid], binApexRts[mid]));
            }
        }
        /// <summary>
        /// Legacy progressive-widening bootstrap. Used as fallback when the center-outward
        /// bootstrap cannot produce enough anchors (e.g. very small libraries, degenerate
        /// RT distributions). Tries windows ±1.0, ±1.5, ±2.0, ±3.0, ±5.0 min centered
        /// on raw library RT values (assumes slope ≈ 1).
        /// </summary>
        private (RtCalibrationModel Model, IRtCalibrationModel DetailedModel,
                 List<CalibrationIterationLog> Log, List<DiaSearchResult> Results)
            LegacyProgressiveBootstrap(
                IList<LibraryPrecursorInput> precursors,
                DiaScanIndex scanIndex,
                DiaSearchParameters baseParameters,
                DiaExtractionOrchestrator orchestrator,
                List<CalibrationIterationLog> log)
        {
            double[] fallbackWindows = { 1.0, 1.5, 2.0, 3.0, 5.0 };
            DiaLibraryQueryGenerator.GenerationResult genResult = default;
            List<DiaSearchResult> bootstrapResults = null;
            bool gotEnough = false;
            var sw = Stopwatch.StartNew();

            foreach (double tryWindow in fallbackWindows)
            {
                Report($"  [Calibration] Legacy bootstrap: trying ±{tryWindow:F1} min window...");
                var bParams = CloneWithRtTolerance(baseParameters, (float)tryWindow);
                genResult = DiaLibraryQueryGenerator.Generate(precursors, scanIndex, bParams);
                var extraction = orchestrator.ExtractAll(genResult.Queries);
                bootstrapResults = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    precursors, genResult, extraction, baseParameters, scanIndex);

                var testAnchors = SelectAnchors(bootstrapResults, precursors, genResult.PrecursorGroups,
                    minApexScore: 0f, InitialTopK, requireMinFragDetRate: false);

                Report($"  [Calibration] Legacy bootstrap ±{tryWindow:F1} min: {testAnchors.Count} anchors");

                if (testAnchors.Count >= MinAnchorCount)
                {
                    gotEnough = true;
                    break;
                }
            }
            sw.Stop();

            if (!gotEnough || bootstrapResults == null)
            {
                log.Add(CreateLogEntry(0, 0, null, baseParameters.RtToleranceMinutes,
                    sw.Elapsed, TimeSpan.Zero, TimeSpan.Zero, "LegacyBootstrapFailed"));
                return (null, null, log, bootstrapResults ?? new List<DiaSearchResult>());
            }

            var anchors = SelectAnchors(bootstrapResults, precursors, genResult.PrecursorGroups,
                minApexScore: 0f, InitialTopK, requireMinFragDetRate: false);

            var model = FitBootstrapModel(anchors, out string modelLabel);
            if (model == null)
            {
                log.Add(CreateLogEntry(0, anchors.Count, null, baseParameters.RtToleranceMinutes,
                    sw.Elapsed, TimeSpan.Zero, TimeSpan.Zero, "LegacyBootstrapFitFailed"));
                return (null, null, log, bootstrapResults);
            }

            double hw = Math.Max(model.SigmaMinutes * SigmaMultiplier, MinWindowHalfWidthMinutes);
            log.Add(CreateLogEntry(0, anchors.Count, model, hw,
                sw.Elapsed, TimeSpan.Zero, TimeSpan.Zero, "Legacy_" + modelLabel));
            Report($"  [Calibration] Legacy bootstrap: slope={model.Slope:F4}, σ={model.SigmaMinutes:F3}, R²={model.RSquared:F4}");

            return RunRefinementIterations(
                precursors, scanIndex, baseParameters, orchestrator,
                model, bootstrapResults, model.SigmaMinutes, model.Slope,
                log, startIteration: 1);
        }

        /// <summary>
        /// Fits the bootstrap linear model from anchors, with RANSAC fallback.
        /// </summary>
        /// <summary>
        /// Fits the bootstrap linear model from anchors using weighted RANSAC.
        /// QualityScore is used as a weight — anchors with lower scores contribute
        /// less by replicating points proportional to their relative weight.
        /// This ensures Phase C down-weighting (×0.3) has real effect on the fit.
        /// </summary>
        private IRtCalibrationModel FitBootstrapModel(
            List<(double LibraryRt, double ObservedRt, double QualityScore)> anchors,
            out string modelLabel)
        {
            modelLabel = "Linear";
            if (anchors.Count == 0) return null;

            // Normalize weights and replicate points so weighted anchors
            // have proportional influence on the RANSAC fit.
            double maxWeight = anchors.Max(a => a.QualityScore);
            if (maxWeight <= 0) maxWeight = 1.0;

            // Scale weights to [0.3, 1.0] range so even down-weighted Phase C
            // anchors contribute at least once. Full-weight anchors replicate 3×.
            const int MaxReplications = 3;
            var libraryRtsList = new List<double>(anchors.Count * MaxReplications);
            var observedRtsList = new List<double>(anchors.Count * MaxReplications);

            foreach (var (libRt, obsRt, quality) in anchors)
            {
                double normalizedWeight = quality / maxWeight;
                // Map [0, 1] → [1, MaxReplications] replications
                int reps = Math.Max(1, (int)Math.Round(normalizedWeight * MaxReplications));
                for (int r = 0; r < reps; r++)
                {
                    libraryRtsList.Add(libRt);
                    observedRtsList.Add(obsRt);
                }
            }

            var libraryRts = libraryRtsList.ToArray();
            var observedRts = observedRtsList.ToArray();

            // With only the center iRT band (narrow range), slope and intercept are
            // highly correlated — RANSAC can fit many (slope, intercept) pairs equally
            // well, consistently returning a wrong slope. Force intercept through the
            // known gradient origin to break the degeneracy.
            // For Prosit CiRT: iRT=0 elutes at ~RT=1 min (gradient start).
            // We estimate the forced intercept from the data rather than hardcoding it.
            double forcedIntercept = EstimateIntercept(libraryRtsList, observedRtsList);
            var forcedModel = FitForcedInterceptLinear(libraryRts, observedRts, forcedIntercept);
            if (forcedModel != null && forcedModel.IsReliable)
            {
                Report($"  [Calibration] ForcedIntercept fit: intercept={forcedIntercept:F3} slope={forcedModel.Slope:F4} σ={forcedModel.SigmaMinutes:F3}");
                return new LinearRtModelWrapper(forcedModel);
            }
            Report($"  [Calibration] ForcedIntercept fit failed (intercept={forcedIntercept:F3}) — falling back to RANSAC");
            var linearModel = FitLinearWithRansac(libraryRts, observedRts, inlierThreshold: 2.0);
            if (linearModel == null || !linearModel.IsReliable)
                linearModel = FitLinearWithRansac(libraryRts, observedRts, inlierThreshold: 3.0);
            if (linearModel == null || !linearModel.IsReliable)
                return null;
            return new LinearRtModelWrapper(linearModel);
        }
        /// <summary>
        /// Estimates the RT-axis intercept by finding the minimum plausible RT
        /// for the lowest-iRT anchors. Uses the 5th percentile of observed RT
        /// values for anchors in the bottom 10% of the iRT range.
        /// Falls back to 1.0 min (typical gradient start) if too few low-iRT anchors.
        /// </summary>
        private static double EstimateIntercept(
            List<double> libraryRts, List<double> observedRts)
        {
            double irtMin = double.MaxValue, irtMax = double.MinValue;
            for (int i = 0; i < libraryRts.Count; i++)
            {
                if (libraryRts[i] < irtMin) irtMin = libraryRts[i];
                if (libraryRts[i] > irtMax) irtMax = libraryRts[i];
            }

            double irtThreshold = irtMin + (irtMax - irtMin) * 0.10;
            var lowIrtObservedRts = new List<double>();
            for (int i = 0; i < libraryRts.Count; i++)
                if (libraryRts[i] <= irtThreshold)
                    lowIrtObservedRts.Add(observedRts[i]);

            if (lowIrtObservedRts.Count < 5)
                return 1.0; // fallback: typical gradient start

            lowIrtObservedRts.Sort();
            return lowIrtObservedRts[lowIrtObservedRts.Count / 20]; // 5th percentile
        }

        /// <summary>
        /// Fits RT = slope * iRT + intercept with the intercept fixed.
        /// Solves for slope only via weighted least squares.
        /// </summary>
        private static RtCalibrationModel FitForcedInterceptLinear(
            double[] libraryRts, double[] observedRts, double intercept)
        {
            if (libraryRts.Length < 5) return null;

            // Shift observed RTs by intercept, then solve slope = sum(x*(y-b)) / sum(x*x)
            double sumXY = 0, sumXX = 0;
            for (int i = 0; i < libraryRts.Count(); i++)
            {
                double x = libraryRts[i];
                double y = observedRts[i] - intercept;
                sumXY += x * y;
                sumXX += x * x;
            }

            if (sumXX < 1e-10) return null;
            double slope = sumXY / sumXX;
            if (slope <= 0) return null;

            // Compute sigma from residuals
            double ssRes = 0;
            for (int i = 0; i < libraryRts.Length; i++)
            {
                double pred = slope * libraryRts[i] + intercept;
                double res = observedRts[i] - pred;
                ssRes += res * res;
            }
            double sigma = Math.Sqrt(ssRes / Math.Max(libraryRts.Length - 2, 1));

            double meanObs = 0;
            for (int i = 0; i < observedRts.Length; i++) meanObs += observedRts[i];
            meanObs /= observedRts.Length;
            double ssTot = 0;
            for (int i = 0; i < observedRts.Length; i++)
            {
                double d = observedRts[i] - meanObs;
                ssTot += d * d;
            }
            double rSquared = ssTot > 0 ? Math.Max(0, 1.0 - ssRes / ssTot) : 0;

            return new RtCalibrationModel(
                slope: slope,
                intercept: intercept,
                sigmaMinutes: Math.Max(sigma, 0.1),
                rSquared: rSquared,
                anchorCount: libraryRts.Length);
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
        /// Selects calibration anchors from extraction results.
        /// Uses iRT (preferred) then RetentionTime as the library RT coordinate.
        /// Anchors are targets with ObservedApexRt and ApexScore >= minApexScore.
        /// </summary>
        public static List<(double LibraryRt, double ObservedRt, double QualityScore)> SelectAnchors(
            List<DiaSearchResult> results,
            IList<LibraryPrecursorInput> precursors,
            DiaLibraryQueryGenerator.PrecursorQueryGroup[] groups,
            float minApexScore,
            int maxAnchors,
            bool requireMinFragDetRate = true,
            Action<string> logger = null)
        {
            var precursorLookup = new Dictionary<(string, int), LibraryPrecursorInput>(
                precursors?.Count ?? 0);
            if (precursors != null)
                foreach (var p in precursors)
                {
                    if (p.IsDecoy) continue;
                    var key = (p.Sequence, p.ChargeState);
                    if (!precursorLookup.ContainsKey(key))
                        precursorLookup[key] = p;
                }

            var candidates = new List<(double LibraryRt, double ObservedRt, double QualityScore)>();

            // Margin to detect boundary-clipped apexes.
            // If ObservedApexRt is within this distance of the extraction window edge,
            // the true apex was likely outside the window and the observed value is biased.
            // One scan cycle ≈ 0.05 min is a safe minimum; use 10% of the window half-width
            // to scale with window size.
            const float BoundaryMarginFraction = 0.10f;
            const float BoundaryMarginMinMinutes = 0.04f;

            foreach (var r in results)
            {
                if (r.IsDecoy) continue;
                if (float.IsNaN(r.ObservedApexRt)) continue;
                if (minApexScore > 0 && (float.IsNaN(r.ApexScore) || r.ApexScore < minApexScore))
                    continue;

                float fragDetRate = r.FragmentDetectionRate;
                if (requireMinFragDetRate && fragDetRate < 0.5f) continue;

                // Reject apex at extraction window boundary — ObservedApexRt is biased.
                float windowHalfWidth = (r.RtWindowEnd - r.RtWindowStart) * 0.5f;
                float margin = Math.Max(windowHalfWidth * BoundaryMarginFraction, BoundaryMarginMinMinutes);
                if (r.ObservedApexRt <= r.RtWindowStart + margin) continue;
                if (r.ObservedApexRt >= r.RtWindowEnd - margin) continue;

                double libRt;
                if (precursorLookup.TryGetValue((r.Sequence, r.ChargeState), out var precursor))
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

                double quality = (minApexScore > 0 ? r.ApexScore : 1f) * fragDetRate * massAccuracyTerm;
                candidates.Add((libRt, r.ObservedApexRt, quality));
            }

            logger?.Invoke($"  [Calibration] SelectAnchors: {candidates.Count} candidates (minApexScore={minApexScore:F2})");

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
                    rtMin = Math.Max((float)(predictedRt - windowHalfWidthMinutes), globalRtMin);
                    rtMax = Math.Min((float)(predictedRt + windowHalfWidthMinutes), globalRtMax);
                }
                else if (p.RetentionTime.HasValue)
                {
                    float predictedRt = (float)model.ToMinutes(p.RetentionTime.Value);
                    rtMin = Math.Max((float)(predictedRt - windowHalfWidthMinutes), globalRtMin);
                    rtMax = Math.Min((float)(predictedRt + windowHalfWidthMinutes), globalRtMax);
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