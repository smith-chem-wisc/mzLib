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
        public double NonLinearSigmaImprovementThreshold { get; set; } = 0.05;

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

            var phaseAMedians = new List<(double Irt, double ApexRt)>();
            var phaseBLoMedians = new List<(double Irt, double ApexRt)>();
            var phaseBHiMedians = new List<(double Irt, double ApexRt)>();
            var phaseCLoMedians = new List<(double Irt, double ApexRt)>();
            var phaseCHiMedians = new List<(double Irt, double ApexRt)>();

            // ── Compute TIC seed model ─────────────────────────────────────────────
            // Maps library iRT extremes to the RT boundaries where 10%/90% of the
            // cumulative MS2 TIC has elapsed. This gives a correct prior slope even
            // when the iRT scale is nonlinear or offset (e.g. Prosit CiRT where
            // iRT=0 elutes at ~16 min, not at run start).
            double irtLow = allTargetsSorted[0].LibRt;
            double irtHigh = allTargetsSorted[allTargetsSorted.Count - 1].LibRt;

            var (ticLowRt, ticHighRt) = EstimateTicBoundaries(scanIndex,
                lowFraction: 0.05, highFraction: 0.95);

            Report($"  [Calibration] TIC boundaries: 10%→{ticLowRt:F2} min, 90%→{ticHighRt:F2} min");
            Report($"  [Calibration] Library iRT range: [{irtLow:F1}, {irtHigh:F1}]");

            // ── TIC diagnostic: report cumulative fractions at known anchor RTs ───
            {
                var ticProfile = scanIndex.GetMs2TicProfile();
                double totalTic = 0;
                for (int i = 0; i < ticProfile.Length; i++)
                    totalTic += ticProfile[i].Intensity;

                if (totalTic > 0)
                {
                    double cumAt11 = 0, cumAt41 = 0;
                    for (int i = 0; i < ticProfile.Length; i++)
                    {
                        if (ticProfile[i].Rt <= 11.0f) cumAt11 += ticProfile[i].Intensity;
                        if (ticProfile[i].Rt <= 41.0f) cumAt41 += ticProfile[i].Intensity;
                    }
                    Report($"  [Calibration] TIC diagnostic: " +
                           $"cumulative at RT=11 min = {cumAt11 / totalTic * 100:F1}%,  " +
                           $"cumulative at RT=41 min = {cumAt41 / totalTic * 100:F1}%");
                    Report($"  [Calibration] Suggested fractions: " +
                           $"lowFraction={cumAt11 / totalTic:F3}, highFraction={cumAt41 / totalTic:F3}");
                }
            }


            double seedSlope = (ticHighRt - ticLowRt) / (irtHigh - irtLow);
            double seedIntercept = ticLowRt - seedSlope * irtLow;
            double seedSigma = Math.Max((ticHighRt - ticLowRt) * 0.08, 1.0); // ~8% of span, ≥1 min

            IRtCalibrationModel ticSeedModel = new LinearRtModelWrapper(new RtCalibrationModel(
                slope: seedSlope,
                intercept: seedIntercept,
                sigmaMinutes: seedSigma,
                rSquared: 0.0,
                anchorCount: 2));

            Report($"  [Calibration] TIC seed model: slope={seedSlope:F4}, intercept={seedIntercept:F3}, " +
                   $"σ={seedSigma:F3} (iRT {irtLow:F1}→RT {ticLowRt:F2}, iRT {irtHigh:F1}→RT {ticHighRt:F2})");

            // ── Phase A: iterate center band (rank 23–43%) ──────────────────────────
            // The TIC seed model is most accurate in the lower-middle iRT range,
            // where the iRT→RT curve is approximately linear.
            // Using rank 40–60% puts the initial search windows ~5 min too early
            // due to the nonlinear compression at high iRT values.
            int centerLo = (int)(totalTargets * 0.23);
            int centerHi = (int)(totalTargets * 0.43);
            var centerIndices = new List<int>();
            for (int i = centerLo; i < centerHi; i++)
                centerIndices.Add(allTargetsSorted[i].Idx);

            var (centerBand, centerExtractionIrts) = BuildBandPrecursorList(precursors, centerIndices, includeDecoys: true);

            IRtCalibrationModel model0 = null;
            double prevSlopeA = double.NaN;
            List<DiaSearchResult> centerResults = null;

            for (int phaseAIter = 0; phaseAIter < 5; phaseAIter++)
            {
                DiaLibraryQueryGenerator.GenerationResult gen;
                double halfWidth;

                if (phaseAIter == 0)
                {
                    halfWidth = Math.Max(ticSeedModel.SigmaMinutes * 3.0, 2.0);
                    gen = GenerateWithModelWindows(
                        centerBand, centerExtractionIrts, ticSeedModel,
                        cellHalfWidth: halfWidth,
                        scanIndex, baseParameters);
                }
                else
                {
                    halfWidth = Math.Max(model0.SigmaMinutes * 2.0, MinWindowHalfWidthMinutes);
                    gen = GenerateWithModelWindows(
                        centerBand, centerExtractionIrts, model0,
                        cellHalfWidth: halfWidth,
                        scanIndex, baseParameters);
                }

                extractSw.Restart();
                var extract = orchestrator.ExtractAll(gen.Queries);
                centerResults = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    centerBand, gen, extract, baseParameters, scanIndex);
                extractSw.Stop();

                // ── Diagnostic: sample apex RTs after extraction ──────────────────
                if (phaseAIter == 0 && centerResults != null && centerResults.Count > 0)
                {
                    int step = Math.Max(1, centerResults.Count / 5);
                    Report($"  [Calibration] Phase A iter 0 apex sample (observedApexRt vs seedPredicted):");
                    for (int di = 0; di < centerResults.Count; di += step)
                    {
                        var r = centerResults[di];
                        double seedPred = ticSeedModel.ToMinutes(r.LibraryRetentionTime ?? 0);
                        Report($"    iRT={r.LibraryRetentionTime:F1}  apexRt={r.ObservedApexRt:F2}  seedPred={seedPred:F2}  isDecoy={r.IsDecoy}");
                    }
                }

                if (phaseAIter == 0 && centerResults != null)
                {
                    var gtPoints = new (double Irt, double TrueRt)[]
                    {
                        (20.1, 19.5), (27.6, 21.5), (33.7, 22.5), (40.6, 23.5),
                        (50.6, 24.5), (61.9, 26.5), (75.2, 28.5), (77.9, 29.5),
                        (81.0, 30.5), (86.1, 31.5), (91.3, 32.5)
                    };

                    Report($"  [GT PROBE] Ground truth center band diagnostics (Phase A iter 0):");
                    foreach (var (irt, trueRt) in gtPoints)
                    {
                        var match = centerResults
                            .Where(r => !r.IsDecoy && r.LibraryRetentionTime.HasValue)
                            .OrderBy(r => Math.Abs(r.LibraryRetentionTime.Value - irt))
                            .FirstOrDefault();

                        if (match == null) continue;

                        string candidates = "none";
                        if (match.DetectedPeakGroup.HasValue && match.DetectedPeakGroup.Value.CandidateApexRts != null)
                            candidates = string.Join(", ", match.DetectedPeakGroup.Value.CandidateApexRts.Select(r => r.ToString("F2")));

                        Report($"  [GT PROBE] iRT={irt:F1} trueRt={trueRt:F1} " +
                               $"seq={match.Sequence} matchIrt={match.LibraryRetentionTime:F1} " +
                               $"apexRt={match.ObservedApexRt:F2} apexScore={match.ApexScore:F3} " +
                               $"window={match.RtWindowStart:F2}-{match.RtWindowEnd:F2} " +
                               $"candidates=[{candidates}]");
                    }
                }

                // ── Diagnostic: candidate count distribution ──────────────────────
                if (phaseAIter == 0 && centerResults != null && centerResults.Count > 0)
                {
                    
                    int singleCandidate = 0, multiCandidate = 0;
                    foreach (var r in centerResults)
                    {
                        if (r.DetectedPeakGroup.HasValue)
                        {
                            int nc = r.DetectedPeakGroup.Value.CandidateCount;
                            if (nc <= 1) singleCandidate++;
                            else multiCandidate++;
                        }
                    }
                    int total = singleCandidate + multiCandidate;
                    Report($"  [Calibration] Phase A iter 0 candidate counts: {singleCandidate} single, {multiCandidate} multi ({(total > 0 ? 100.0 * multiCandidate / total : 0):F0}% multi)");
                }
                fitSw.Restart();
                var candidate = ComputeGridModel(centerResults, centerBand,
    accumulator: phaseAMedians,
    rtPercentile: phaseAIter == 0 ? 0.85 : 0.50);
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

                if (!double.IsNaN(prevSlopeA) && Math.Abs(candidate.Slope - prevSlopeA) / prevSlopeA < 0.02)
                {
                    Report($"  [Calibration] Phase A converged at iter {phaseAIter}.");
                    break;
                }
                prevSlopeA = candidate.Slope;
            }

            Report($"  [Calibration] Phase A final: slope={model0.Slope:F4}, intercept={model0.Intercept:F3}, σ={model0.SigmaMinutes:F3}");

            // ── Phase B: mid bands (rank 20–40%, 60–80%) → model_1 ────────────────
            IRtCalibrationModel modelBLo = null;
            IRtCalibrationModel modelBHi = null;

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

            var (midBand, midExtractionIrts) = BuildBandPrecursorList(precursors, midIndices, includeDecoys: true);
            var (midBandLo, _) = BuildBandPrecursorList(precursors, midIndicesLo, includeDecoys: false);
            var (midBandHi, _) = BuildBandPrecursorList(precursors, midIndicesHi, includeDecoys: false);
            var midLoKeys = new HashSet<(string, int)>(midBandLo.Select(p => (p.Sequence, p.ChargeState)));

            IRtCalibrationModel model1 = model0;
            double prevSlopeB = double.NaN;
            List<DiaSearchResult> midResults = null;

            for (int phaseBIter = 0; phaseBIter < 5; phaseBIter++)
            {
                phaseBLoMedians.Clear();
                phaseBHiMedians.Clear();

                double halfWidth = phaseBIter == 0
                    ? Math.Max(model1.SigmaMinutes * 3.0, 2.0)
                    : Math.Max(model1.SigmaMinutes * 2.0, MinWindowHalfWidthMinutes);

                var gen = GenerateWithModelWindows(
                    midBand, midExtractionIrts, model1,
                    cellHalfWidth: halfWidth,
                    scanIndex, baseParameters);

                extractSw.Restart();
                var extract = orchestrator.ExtractAll(gen.Queries);
                midResults = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    midBand, gen, extract, baseParameters, scanIndex);
                extractSw.Stop();

                var midResultsLo = midResults.Where(r => midLoKeys.Contains((r.Sequence, r.ChargeState))).ToList();
                var midResultsHi = midResults.Where(r => !midLoKeys.Contains((r.Sequence, r.ChargeState))).ToList();

                fitSw.Restart();
                var candidate = ComputeGridModelSplit(
                    midResultsLo, midBandLo, midResultsHi, midBandHi,
                    accumulator: null);
                var candidateLo = ComputeGridModel(midResultsLo, midBandLo, accumulator: phaseBLoMedians);
                var candidateHi = ComputeGridModel(midResultsHi, midBandHi, accumulator: phaseBHiMedians);
                if (candidateLo != null) modelBLo = candidateLo;
                if (candidateHi != null) modelBHi = candidateHi;
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
            IRtCalibrationModel modelCLo = modelBLo ?? model1;
            IRtCalibrationModel modelCHi = modelBHi ?? model1;

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

            var (outerBand, outerExtractionIrts) = BuildBandPrecursorList(precursors, outerIndices, includeDecoys: true);
            var (outerBandLo, _) = BuildBandPrecursorList(precursors, outerIndicesLo, includeDecoys: false);
            var (outerBandHi, _) = BuildBandPrecursorList(precursors, outerIndicesHi, includeDecoys: false);
            var outerLoKeys = new HashSet<(string, int)>(outerBandLo.Select(p => (p.Sequence, p.ChargeState)));

            IRtCalibrationModel modelC = model1;
            double prevSlopeC = double.NaN;
            List<DiaSearchResult> outerResults = null;
            List<DiaSearchResult> outerResultsLo = null;
            List<DiaSearchResult> outerResultsHi = null;
            for (int phaseCIter = 0; phaseCIter < 5; phaseCIter++)
            {
                phaseCLoMedians.Clear();
                phaseCHiMedians.Clear();

                double halfWidthLo = Math.Max(modelCLo.SigmaMinutes * 2.0, MinWindowHalfWidthMinutes);
                double halfWidthHi = Math.Max(modelCHi.SigmaMinutes * 2.0, MinWindowHalfWidthMinutes);

                var (outerBandLoWithDecoys, outerExtractionIrtsLo) = BuildBandPrecursorList(precursors, outerIndicesLo, includeDecoys: true);
                var (outerBandHiWithDecoys, outerExtractionIrtsHi) = BuildBandPrecursorList(precursors, outerIndicesHi, includeDecoys: true);

                var genLo = GenerateWithModelWindows(outerBandLoWithDecoys, outerExtractionIrtsLo, modelCLo,
                    cellHalfWidth: halfWidthLo, scanIndex, baseParameters);
                var genHi = GenerateWithModelWindows(outerBandHiWithDecoys, outerExtractionIrtsHi, modelCHi,
                    cellHalfWidth: halfWidthHi, scanIndex, baseParameters);

                extractSw.Restart();
                var extractLo = orchestrator.ExtractAll(genLo.Queries);
                var extractHi = orchestrator.ExtractAll(genHi.Queries);
                outerResultsLo = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    outerBandLoWithDecoys, genLo, extractLo, baseParameters, scanIndex);
                outerResultsHi = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    outerBandHiWithDecoys, genHi, extractHi, baseParameters, scanIndex);
                extractSw.Stop();

                var candidateLo = ComputeGridModel(outerResultsLo, outerBandLo, accumulator: phaseCLoMedians);
                var candidateHi = ComputeGridModel(outerResultsHi, outerBandHi, accumulator: phaseCHiMedians);

                var candidate = ComputeGridModelSplit(
                    outerResultsLo, outerBandLo, outerResultsHi, outerBandHi,
                    accumulator: null);

                if (candidate == null)
                {
                    Report($"  [Calibration] Phase C iter {phaseCIter} grid failed — keeping previous model.");
                    break;
                }

                if (candidateLo != null) modelCLo = candidateLo;
                if (candidateHi != null) modelCHi = candidateHi;

                Report($"  [Calibration] Phase C iter {phaseCIter}: slope={candidate.Slope:F4}, intercept={candidate.Intercept:F3}, σ={candidate.SigmaMinutes:F3}, R²={candidate.RSquared:F4} (lo±{halfWidthLo:F2} hi±{halfWidthHi:F2} min)");

                if (!double.IsNaN(prevSlopeC) && Math.Abs(candidate.Slope - prevSlopeC) / prevSlopeC < 0.02)
                {
                    Report($"  [Calibration] Phase C converged at iter {phaseCIter}: lo={phaseCLoMedians.Count}, hi={phaseCHiMedians.Count}");
                    break;
                }
                prevSlopeC = candidate.Slope;
            }

            // ── Select final bootstrap model and cap sigma ────────────────────────
            const double BootstrapSigmaCap = 0.3;
            IRtCalibrationModel bootstrapFinal = modelC;

            var allPhaseMedians = phaseAMedians
                .Concat(phaseBLoMedians)
                .Concat(phaseBHiMedians)
                .Concat(phaseCLoMedians)
                .Concat(phaseCHiMedians)
                .OrderBy(m => m.Irt)
                .ToList();

            Report($"  [Calibration] Bootstrap model selection: A={phaseAMedians.Count}, BLo={phaseBLoMedians.Count}, BHi={phaseBHiMedians.Count}, CLo={phaseCLoMedians.Count}, CHi={phaseCHiMedians.Count}, total={allPhaseMedians.Count}, modelC.σ={modelC.SigmaMinutes:F4}, EnableNonLinear={EnableNonLinearModelSelection}");

            if (EnableNonLinearModelSelection && allPhaseMedians.Count >= 50)
            {
                var libRts = allPhaseMedians.Select(m => m.Irt).ToArray();
                var apexRts = allPhaseMedians.Select(m => m.ApexRt).ToArray();

                var lowessModel = LowessRtModel.Fit(libRts, apexRts, bandwidth: 0.25, enforceMonotonic: true);
                Report($"  [Calibration] LOWESS fit: null={lowessModel == null}, reliable={lowessModel?.IsReliable}, σ={lowessModel?.SigmaMinutes:F4}, threshold={modelC.SigmaMinutes:F4}"); 
                Report($"  [Calibration] LOWESS predictions: iRT=0→{lowessModel.ToMinutes(0):F2}, iRT=34→{lowessModel.ToMinutes(34):F2}, iRT=59→{lowessModel.ToMinutes(59):F2}, iRT=84→{lowessModel.ToMinutes(84):F2}, iRT=110→{lowessModel.ToMinutes(110):F2}, iRT=170→{lowessModel.ToMinutes(170):F2}");
                if (lowessModel != null && lowessModel.IsReliable)
                {
                    Report($"  [Calibration] LOWESS bootstrap model selected: σ={lowessModel.SigmaMinutes:F3}, R²={lowessModel.RSquared:F4}");
                    // Use BootstrapSigmaCap for window generation — LOWESS σ is fit residual on
                    // clean bin medians, not real anchor scatter
                    bootstrapFinal = new LowessSigmaOverrideModel(lowessModel, BootstrapSigmaCap);
                }
                else
                {
                    // Fall back to piecewise linear
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
            }

            // Cap sigma for window generation but preserve the nonlinear curve shape
            if (bootstrapFinal.SigmaMinutes > BootstrapSigmaCap)
            {
                // IMPORTANT: if we have a nonlinear model, keep it — only override sigma
                // by wrapping it in a SigmaOverrideModel rather than collapsing to linear
                if (bootstrapFinal is LowessRtModel lowess)
                {
                    currentModel = new LowessSigmaOverrideModel(lowess, BootstrapSigmaCap);
                    Report($"  [Calibration] Bootstrap complete (LOWESS, σ capped): σ={bootstrapFinal.SigmaMinutes:F3}→{BootstrapSigmaCap:F3}, R²={bootstrapFinal.RSquared:F4}");
                }
                else
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
            }
            else
            {
                currentModel = bootstrapFinal;
                Report($"  [Calibration] Bootstrap complete: slope={currentModel.Slope:F4}, intercept={currentModel.Intercept:F3}, σ={currentModel.SigmaMinutes:F3}, R²={currentModel.RSquared:F4}");
            }

            // Combine all results for refinement iterations
            currentResults = new List<DiaSearchResult>(centerResults.Count + midResults.Count + outerResultsLo.Count + outerResultsHi.Count);
            currentResults.AddRange(centerResults);
            currentResults.AddRange(midResults);
            currentResults.AddRange(outerResultsLo);
            currentResults.AddRange(outerResultsHi);

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
        /// <summary>
        /// Finds the RT at which cumulative MS2 TIC crosses lowFraction and highFraction.
        /// Uses scan peak sums as TIC proxy. Falls back to trimmed file RT range if sparse.
        /// </summary>
        private (double LowRt, double HighRt) EstimateTicBoundaries(
            DiaScanIndex scanIndex,
            double lowFraction = 0.10,
            double highFraction = 0.90)
        {
            var profile = scanIndex.GetMs2TicProfile(); // (float Rt, float Intensity)[]

            // Fallback: use trimmed file RT range
            double fallbackLow = scanIndex.GetGlobalRtMin() + 2.0;
            double fallbackHigh = scanIndex.GetGlobalRtMax() - 2.0;

            if (profile == null || profile.Length < 10)
                return (fallbackLow, fallbackHigh);

            double total = 0;
            for (int i = 0; i < profile.Length; i++)
                total += profile[i].Intensity;

            if (total <= 0) return (fallbackLow, fallbackHigh);

            double cumulative = 0;
            double lowRt = profile[0].Rt;
            double highRt = profile[profile.Length - 1].Rt;
            bool foundLow = false;

            for (int i = 0; i < profile.Length; i++)
            {
                cumulative += profile[i].Intensity;
                double fraction = cumulative / total;

                if (!foundLow && fraction >= lowFraction)
                {
                    lowRt = profile[i].Rt;
                    foundLow = true;
                }
                if (fraction >= highFraction)
                {
                    highRt = profile[i].Rt;
                    break;
                }
            }

            // Sanity check — must be at least 5 min apart
            if (highRt - lowRt < 5.0)
                return (fallbackLow, fallbackHigh);

            return (lowRt, highRt);
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

                fitSw.Stop();

                double newSigma = bestModel.SigmaMinutes;
                double newSlope = bestModel.Slope;

                Console.WriteLine($"  [Calibration] Iteration {iteration}: slope={newSlope:F4}, σ={newSigma:F3}, R²={bestModel.RSquared:F4}, model={modelTypeLabel}");

                double divergenceThreshold = iteration <= 1 ? 1.20 : 1.10;
                if (newSigma > previousSigma * divergenceThreshold)
                {
                    Console.WriteLine($"  [Calibration] Iteration {iteration}: σ diverged ({newSigma:F3} > {previousSigma * 1.1:F3}), reverting.");
                    log.Add(CreateLogEntry(iteration, anchors.Count, bestModel,
                        Math.Min(Math.Max(newSigma * SigmaMultiplier, MinRefinementWindowMinutes), maxWindowHalfWidth),
                        extractSw.Elapsed, fitSw.Elapsed, selectSw.Elapsed, modelTypeLabel + "_Diverged_Reverted"));
                    break;
                }

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
                double nextRawWindow = Math.Max(newSigma * SigmaMultiplier, MinRefinementWindowMinutes);
                if (nextRawWindow < maxWindowHalfWidth)
                    maxWindowHalfWidth = nextRawWindow;
            }

            if (currentModel == null)
                return (null, null, log, currentResults);

            return (currentModel.ToRtCalibrationModel(), currentModel, log, currentResults);
        }

        // ── Phase 23 Bootstrap Helpers ─────────────────────────────────────

        private (double Min, double Max) EstimatePeptideElutionWindow(DiaScanIndex scanIndex)
        {
            double rtMin = scanIndex.GetGlobalRtMin();
            double rtMax = scanIndex.GetGlobalRtMax();
            double span = rtMax - rtMin;

            const double TrimMinutes = 2.0;

            if (span > TrimMinutes * 4)
            {
                rtMin += TrimMinutes;
                rtMax -= TrimMinutes;
            }

            return (rtMin, rtMax);
        }

        /// <summary>
        /// Builds a precursor list for a single iRT band with paired extraction iRTs.
        /// Targets and their paired decoys are interleaved: (target0, decoy0, target1, decoy1, ...).
        /// extractionIrts[i] contains the iRT to use for window placement of precursors[i].
        /// For targets: their own iRT. For decoys: the paired target's iRT.
        /// This ensures each decoy competes in the same local signal environment as its target.
        /// </summary>
        private static (List<LibraryPrecursorInput> precursors, double?[] extractionIrts)
            BuildBandPrecursorList(
                IList<LibraryPrecursorInput> all,
                IList<int> targetIndices,
                bool includeDecoys)
        {
            int firstDecoyIdx = -1;
            if (includeDecoys)
            {
                for (int i = 0; i < all.Count; i++)
                {
                    if (all[i].IsDecoy) { firstDecoyIdx = i; break; }
                }
            }

            var result = new List<LibraryPrecursorInput>(
                targetIndices.Count * (includeDecoys ? 2 : 1));
            var extractionIrts = new List<double?>();

            foreach (int idx in targetIndices)
            {
                var target = all[idx];
                double? targetIrt = target.IrtValue ?? target.RetentionTime;

                result.Add(target);
                extractionIrts.Add(targetIrt);  // target uses its own iRT

                if (!includeDecoys || firstDecoyIdx < 0)
                    continue;

                int decoyIdx = firstDecoyIdx + idx;
                if (decoyIdx < all.Count && all[decoyIdx].IsDecoy)
                {
                    result.Add(all[decoyIdx]);
                    extractionIrts.Add(targetIrt);  // decoy uses paired target's iRT
                }
            }

            return (result, extractionIrts.ToArray());
        }

        /// <summary>
        /// Generates queries using a calibration model and per-precursor extraction iRTs.
        /// extractionIrts[i] is the iRT to use for window placement of bandPrecursors[i].
        /// For targets this is their own iRT; for decoys it is the paired target's iRT.
        /// </summary>
        private static DiaLibraryQueryGenerator.GenerationResult GenerateWithModelWindows(
            IList<LibraryPrecursorInput> bandPrecursors,
            double?[] extractionIrts,
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

                // Use the pre-computed extraction iRT (target's own iRT for targets,
                // paired target's iRT for decoys).
                double? irt = (extractionIrts != null && i < extractionIrts.Length)
                    ? extractionIrts[i]
                    : (p.IrtValue ?? p.RetentionTime);

                float rtMin, rtMax;
                if (irt.HasValue && !double.IsNaN(irt.Value))
                {
                    double predictedRt = model.ToMinutes(irt.Value);
                    rtMin = (float)Math.Max(predictedRt - cellHalfWidth, globalRtMin);
                    rtMax = (float)Math.Min(predictedRt + cellHalfWidth, globalRtMax);
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

        private IRtCalibrationModel ComputeGridModel(
            List<DiaSearchResult> results,
            IList<LibraryPrecursorInput> bandPrecursors,
            int targetsPerBin = 50,
            List<(double Irt, double ApexRt)> accumulator = null,
            double rtPercentile = 0.85)
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

            // Collect candidate RTs with their apex cosine scores
            var candidatesByKey = new Dictionary<(string, int), List<(float Rt, float Score)>>();
            foreach (var r in results)
            {
                if (r.IsDecoy) continue;
                var key = (r.Sequence, r.ChargeState);
                if (!candidatesByKey.TryGetValue(key, out var list))
                {
                    list = new List<(float, float)>();
                    candidatesByKey[key] = list;
                }

                if (r.DetectedPeakGroup.HasValue
                    && r.DetectedPeakGroup.Value.CandidateApexRts != null
                    && r.DetectedPeakGroup.Value.CandidateApexRts.Length > 0)
                {
                    // The selected apex has a known score; other candidates get a discounted score
                    float bestRt = r.DetectedPeakGroup.Value.ApexRt;
                    float bestScore = r.ApexScore;
                    foreach (float rt in r.DetectedPeakGroup.Value.CandidateApexRts)
                    {
                        // Selected apex gets full score; others get 50% — they're unscored candidates
                        float score = Math.Abs(rt - bestRt) < 0.05f ? bestScore : bestScore * 0.5f;
                        list.Add((rt, score));
                    }
                }
                else if (!float.IsNaN(r.ObservedApexRt))
                {
                    list.Add((r.ObservedApexRt, r.ApexScore));
                }
            }

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
                var allBinApexRts = new List<(float Rt, float Score)>();

                for (int i = lo; i < hi; i++)
                {
                    var p = targetPrecursors[i];
                    var aKey = (p.Sequence, p.ChargeState);
                    if (!irtByKey.TryGetValue(aKey, out double irt)) continue;
                    if (!candidatesByKey.TryGetValue(aKey, out var candidates)) continue;
                    if (candidates.Count == 0) continue;

                    binIrts.Add(irt);
                    allBinApexRts.AddRange(candidates);
                }

                if (binIrts.Count < 5 || allBinApexRts.Count < 5) continue;

                binIrts.Sort();
                double medianIrt = binIrts[binIrts.Count / 2];

                // Use 85th percentile of all candidate RTs.
                // Interference peaks cluster early; true peptide peaks elute later.
                // The 85th percentile skips past the dense early interference cluster.
                var sortedRts = allBinApexRts.Select(x => x.Rt).OrderBy(x => x).ToList();
                int p85idx = Math.Min(sortedRts.Count - 1, (int)(sortedRts.Count * rtPercentile));
                double medianRt = sortedRts[p85idx];

                binMedianIrt.Add(medianIrt);
                binMedianRt.Add(medianRt);
                // After collecting allBinApexRts for the first bin only:
                if (b == 0)
                {
                    float minCand = allBinApexRts.Min(x => x.Rt);
                    float maxCand = allBinApexRts.Max(x => x.Rt);
                    float meanCand = allBinApexRts.Average(x => x.Rt);
                    Report($"  [DIAG] Bin 0 candidates: n={allBinApexRts.Count} RT=[{minCand:F2},{maxCand:F2}] mean={meanCand:F2}");
                }
            }

            if (accumulator != null)
            {
                for (int i = 0; i < binMedianIrt.Count; i++)
                    accumulator.Add((binMedianIrt[i], binMedianRt[i]));
            }

            Report($"  [Calibration] Apex grid: {numBins} bins, {binMedianIrt.Count} valid (targetsPerBin={targetsPerBin})");

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

            // For each precursor, collect all candidate apex RTs with their scores
            var candidatesByKey = new Dictionary<(string, int), List<(float Rt, float Score)>>();
            foreach (var r in results)
            {
                if (r.IsDecoy) continue;
                var key = (r.Sequence, r.ChargeState);
                if (!candidatesByKey.TryGetValue(key, out var list))
                {
                    list = new List<(float, float)>();
                    candidatesByKey[key] = list;
                }
                if (r.DetectedPeakGroup.HasValue
                    && r.DetectedPeakGroup.Value.CandidateApexRts != null
                    && r.DetectedPeakGroup.Value.CandidateApexRts.Length > 0)
                {
                    float bestRt = r.DetectedPeakGroup.Value.ApexRt;
                    float bestScore = r.ApexScore;
                    foreach (float rt in r.DetectedPeakGroup.Value.CandidateApexRts)
                    {
                        float score = Math.Abs(rt - bestRt) < 0.05f ? bestScore : bestScore * 0.5f;
                        list.Add((rt, score));
                    }
                }
                else if (!float.IsNaN(r.ObservedApexRt))
                {
                    list.Add((r.ObservedApexRt, r.ApexScore));
                }
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
                var binBestRts = new List<double>();

                for (int i = lo; i < hi; i++)
                {
                    var p = targetPrecursors[i];
                    var aKey = (p.Sequence, p.ChargeState);
                    if (!irtByKey.TryGetValue(aKey, out double irt)) continue;
                    if (!candidatesByKey.TryGetValue(aKey, out var candidates)) continue;
                    if (candidates.Count == 0) continue;

                    // Best-scoring candidate for this precursor
                    float bestRt = candidates.OrderByDescending(x => x.Score).First().Rt;
                    binIrts.Add(irt);
                    binBestRts.Add(bestRt);
                }

                if (binIrts.Count < 5) continue;

                binIrts.Sort();
                binBestRts.Sort();
                int mid = binIrts.Count / 2;
                output.Add((binIrts[mid], binBestRts[mid]));
            }
        }

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

        private IRtCalibrationModel FitBootstrapModel(
            List<(double LibraryRt, double ObservedRt, double QualityScore)> anchors,
            out string modelLabel)
        {
            modelLabel = "Linear";
            if (anchors.Count == 0) return null;

            double maxWeight = anchors.Max(a => a.QualityScore);
            if (maxWeight <= 0) maxWeight = 1.0;

            const int MaxReplications = 3;
            var libraryRtsList = new List<double>(anchors.Count * MaxReplications);
            var observedRtsList = new List<double>(anchors.Count * MaxReplications);

            foreach (var (libRt, obsRt, quality) in anchors)
            {
                double normalizedWeight = quality / maxWeight;
                int reps = Math.Max(1, (int)Math.Round(normalizedWeight * MaxReplications));
                for (int r = 0; r < reps; r++)
                {
                    libraryRtsList.Add(libRt);
                    observedRtsList.Add(obsRt);
                }
            }

            var libraryRts = libraryRtsList.ToArray();
            var observedRts = observedRtsList.ToArray();

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
                return 1.0;

            lowIrtObservedRts.Sort();
            return lowIrtObservedRts[lowIrtObservedRts.Count / 20];
        }

        private static RtCalibrationModel FitForcedInterceptLinear(
            double[] libraryRts, double[] observedRts, double intercept)
        {
            if (libraryRts.Length < 5) return null;

            double sumXY = 0, sumXX = 0;
            for (int i = 0; i < libraryRts.Length; i++)
            {
                double x = libraryRts[i];
                double y = observedRts[i] - intercept;
                sumXY += x * y;
                sumXX += x * x;
            }

            if (sumXX < 1e-10) return null;
            double slope = sumXY / sumXX;
            if (slope <= 0) return null;

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

        

        // ── Anchor Selection ───────────────────────────────────────────────

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

        /// <summary>
        /// Generates queries for the full precursor list using a calibration model and
        /// explicit window half-width. Used for refinement iterations.
        /// Decoys are searched at their paired target's predicted RT so they compete
        /// in the same local signal environment (correct T/D competition).
        /// Assembly uses each precursor's own IrtValue for LibraryRetentionTime,
        /// so RtDeviationMinutes is computed correctly for both targets and decoys.
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

            // Find first decoy index for paired-target window placement
            int firstDecoyIdx = -1;
            for (int i = 0; i < precursors.Count; i++)
            {
                if (precursors[i].IsDecoy) { firstDecoyIdx = i; break; }
            }

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

                // For decoys, use the paired target's iRT for window placement.
                // The paired target is at index (i - firstDecoyIdx) in the targets block.
                double? extractionIrt;
                if (p.IsDecoy && firstDecoyIdx >= 0)
                {
                    int pairedTargetIdx = i - firstDecoyIdx;
                    extractionIrt = pairedTargetIdx >= 0 && pairedTargetIdx < firstDecoyIdx
                        ? precursors[pairedTargetIdx].IrtValue ?? precursors[pairedTargetIdx].RetentionTime
                        : p.IrtValue ?? p.RetentionTime;
                }
                else
                {
                    extractionIrt = p.IrtValue ?? p.RetentionTime;
                }

                float rtMin, rtMax;
                if (extractionIrt.HasValue)
                {
                    float predictedRt = (float)model.ToMinutes(extractionIrt.Value);
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