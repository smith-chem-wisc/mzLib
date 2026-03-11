// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// Phase 15, Prompt 4: High-level calibration pipeline orchestration
// Placement: MassSpectrometry/Dia/Calibration/DiaCalibrationPipeline.cs

using System;
using System.Collections.Generic;
using System.Diagnostics;
using MassSpectrometry.Dia.Calibration;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// High-level orchestration for the complete DIA analysis pipeline with automatic
    /// iterative RT calibration.
    ///
    /// Replaces the manual broad → calibrate → narrow pattern with a single method call
    /// that handles:
    ///   1. Iterative calibration (via <see cref="IterativeRtCalibrator"/>)
    ///   2. Final extraction with RT-adaptive windows (via <see cref="DiaLibraryQueryGenerator.GenerateCalibratedAdaptive"/>)
    ///   3. Result assembly with temporal scoring
    ///   4. RT deviation recalibration for downstream feature extraction / FDR
    ///
    /// This class is stateless — all state lives in the parameters and dependencies.
    /// Thread safety: not designed for concurrent calls (the orchestrator and scan index
    /// are shared mutable state), but each call is self-contained.
    /// </summary>
    public static class DiaCalibrationPipeline
    {
        /// <summary>
        /// Result of the full automatic calibration pipeline.
        /// </summary>
        public readonly struct PipelineResult
        {
            /// <summary>Fully scored results ready for feature extraction and FDR.</summary>
            public readonly List<DiaSearchResult> Results;

            /// <summary>The sealed production calibration model (for query generation, FDR).</summary>
            public readonly RtCalibrationModel Calibration;

            /// <summary>The detailed calibration model with GetLocalSigma() for adaptive windows.</summary>
            public readonly IRtCalibrationModel DetailedCalibration;

            /// <summary>Per-iteration calibration diagnostics.</summary>
            public readonly List<CalibrationIterationLog> CalibrationLog;

            /// <summary>Time spent in the calibration loop.</summary>
            public readonly TimeSpan CalibrationTime;

            /// <summary>Time spent in the final extraction pass.</summary>
            public readonly TimeSpan FinalExtractionTime;

            /// <summary>Time spent recalibrating RT deviations.</summary>
            public readonly TimeSpan RtDeviationRecalibrationTime;

            /// <summary>Total pipeline time (calibration + final extraction + recalibration).</summary>
            public TimeSpan TotalTime => CalibrationTime + FinalExtractionTime + RtDeviationRecalibrationTime;

            public PipelineResult(
                List<DiaSearchResult> results,
                RtCalibrationModel calibration,
                IRtCalibrationModel detailedCalibration,
                List<CalibrationIterationLog> calibrationLog,
                TimeSpan calibrationTime,
                TimeSpan finalExtractionTime,
                TimeSpan rtDeviationRecalibrationTime)
            {
                Results = results;
                Calibration = calibration;
                DetailedCalibration = detailedCalibration;
                CalibrationLog = calibrationLog;
                CalibrationTime = calibrationTime;
                FinalExtractionTime = finalExtractionTime;
                RtDeviationRecalibrationTime = rtDeviationRecalibrationTime;
            }
        }

        /// <summary>
        /// Complete DIA analysis pipeline with automatic iterative RT calibration.
        /// Replaces the manual broad → calibrate → narrow pattern.
        ///
        /// Steps:
        ///   1. Run <see cref="IterativeRtCalibrator.Calibrate"/> to converge on a calibration model
        ///   2. Final extraction with <see cref="DiaLibraryQueryGenerator.GenerateCalibratedAdaptive"/>
        ///      using the converged model and RT-adaptive window widths
        ///   3. Assemble results with temporal scoring
        ///   4. Recalibrate RT deviations (RtDeviationMinutes, RtDeviationSquared) using
        ///      the calibration model so downstream classifiers see calibrated deviations
        ///
        /// Returns fully scored results ready for FDR analysis.
        /// </summary>
        /// <param name="precursors">Target + decoy library precursors.</param>
        /// <param name="scanIndex">The DIA scan index (immutable, shared).</param>
        /// <param name="parameters">
        /// Search parameters. RtToleranceMinutes sets the initial broad window for iteration 0.
        /// </param>
        /// <param name="orchestrator">The extraction orchestrator.</param>
        /// <param name="calibrator">
        /// Calibrator with desired settings. Null = use default IterativeRtCalibrator.
        /// </param>
        /// <returns>Pipeline result with scored results, calibration model, and diagnostics.</returns>
        public static PipelineResult RunWithAutomaticCalibration(
            IList<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            DiaSearchParameters parameters,
            DiaExtractionOrchestrator orchestrator,
            IterativeRtCalibrator calibrator = null,
            Action<string> progressReporter = null)
        {
            if (precursors == null) throw new ArgumentNullException(nameof(precursors));
            if (scanIndex == null) throw new ArgumentNullException(nameof(scanIndex));
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));
            if (orchestrator == null) throw new ArgumentNullException(nameof(orchestrator));

            calibrator ??= new IterativeRtCalibrator();

            // Wire progress reporting so calibration messages appear in the GUI/log
            if (progressReporter != null)
                calibrator.ProgressReporter = progressReporter;

            // ── Step 1: Iterative Calibration ──────────────────────────────
            var calibSw = Stopwatch.StartNew();

            var (model, detailedModel, log, calibResults) =
                calibrator.Calibrate(precursors, scanIndex, parameters, orchestrator);

            calibSw.Stop();

            if (model == null)
            {
                // Calibration failed entirely — return whatever results we have
                // from the last iteration (likely the broad-window bootstrap)
                return new PipelineResult(
                    calibResults ?? new List<DiaSearchResult>(),
                    null, null, log,
                    calibSw.Elapsed, TimeSpan.Zero, TimeSpan.Zero);
            }

            // ── Step 2: Final Extraction with Adaptive Windows ─────────────
            // This is a deliberate precision pass using a TIGHTER window than the
            // calibration iterations. The calibration loop uses ±0.506 min (capped by
            // BootstrapSigmaMultiplier=2.0 × σ_bootstrap=0.253) to ensure all real anchors
            // are found. The final re-extraction uses max(SigmaMultiplier×σ_final, MinWindow)
            // = max(4.0×0.077, 0.3) = ±0.31 min — significantly tighter.
            //
            // Tighter final window → fewer decoys → cleaner T/D separation → more IDs at FDR.
            // DO NOT replace this with calibResults directly: that uses the wider calibration
            // window and produces ~3,000 fewer NN IDs at 1% FDR.
            var extractSw = Stopwatch.StartNew();

            var genResult = DiaLibraryQueryGenerator.GenerateCalibratedAdaptive(
                precursors, scanIndex, parameters,
                model, detailedModel,
                sigmaMultiplier: calibrator.SigmaMultiplier,
                minWindowHalfWidthMinutes: calibrator.MinWindowHalfWidthMinutes);

            var extractionResult = orchestrator.ExtractAll(genResult.Queries);

            // ── Step 3: Assemble Results with Temporal Scoring ──────────────
            var results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                precursors, genResult, extractionResult, parameters, scanIndex);

            // ── Step 3b: Chimeric Score ─────────────────────────────────────
            DiaLibraryQueryGenerator.ComputeChimericScores(
                precursors, results, (float)parameters.PpmTolerance);

            extractSw.Stop();

            // ── Step 4: Recalibrate RT Deviations ──────────────────────────
            var rtDevSw = Stopwatch.StartNew();

            RecalibrateRtDeviations(results, precursors, detailedModel);

            rtDevSw.Stop();

            return new PipelineResult(
                results, model, detailedModel, log,
                calibSw.Elapsed, extractSw.Elapsed, rtDevSw.Elapsed);
        }

        /// <summary>
        /// Simplified overload that returns results, calibration, and log as a tuple.
        /// Provided for backward compatibility with the signature specified in the prompt.
        /// </summary>
        public static (List<DiaSearchResult> Results, RtCalibrationModel Calibration,
                        List<CalibrationIterationLog> CalibrationLog)
            RunWithAutomaticCalibrationSimple(
                IList<LibraryPrecursorInput> precursors,
                DiaScanIndex scanIndex,
                DiaSearchParameters parameters,
                DiaExtractionOrchestrator orchestrator,
                IterativeRtCalibrator calibrator = null,
                Action<string> progressReporter = null)
        {
            var result = RunWithAutomaticCalibration(
                precursors, scanIndex, parameters, orchestrator, calibrator, progressReporter);
            return (result.Results, result.Calibration, result.CalibrationLog);
        }

        // ════════════════════════════════════════════════════════════════
        //  Task 3: RT Deviation Recalibration
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Recomputes RtDeviationMinutes and RtDeviationSquared using calibrated predictions
        /// instead of raw library RT.
        ///
        /// Before calibration, RtDeviationMinutes = |ObservedApexRt - LibraryRetentionTime|.
        /// This is problematic because the classifier sees similar deviations for targets
        /// and decoys (both are large and noisy before calibration).
        ///
        /// After calibration:
        ///   predictedRt = calibration.ToMinutes(libraryRt)
        ///   RtDeviationMinutes = ObservedApexRt - predictedRt  (signed, not absolute)
        ///   RtDeviationSquared = RtDeviationMinutes²
        ///
        /// Signed deviation preserves directionality — peptides that elute earlier vs later
        /// than predicted can carry different discriminative information.
        ///
        /// This must be called after result assembly and before feature extraction / FDR.
        /// </summary>
        /// <param name="results">Results to recalibrate. Modified in place.</param>
        /// <param name="precursors">
        /// Precursor inputs, parallel by InputIndex (from PrecursorQueryGroup).
        /// Used to retrieve library RT for each result.
        /// </param>
        /// <param name="calibration">The fitted RT calibration model.</param>
        public static void RecalibrateRtDeviations(
            List<DiaSearchResult> results,
            IList<LibraryPrecursorInput> precursors,
            IRtCalibrationModel calibration)
        {
            if (results == null || calibration == null)
                return;

            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];
                if (float.IsNaN(r.ObservedApexRt))
                    continue;

                double? libRt = r.LibraryRetentionTime;
                if (!libRt.HasValue || double.IsNaN(libRt.Value))
                    continue;

                // Always convert through the calibration model.
                // For iRT libraries: libRt is raw iRT units, ToMinutes applies slope+intercept.
                // For native-RT libraries: libRt is minutes, ToMinutes applies any systematic
                // offset/nonlinearity the model captured. Either way this is the correct
                // predicted RT to compare against ObservedApexRt.
                double predictedRt = calibration.ToMinutes(libRt.Value);

                float deviation = (float)(r.ObservedApexRt - predictedRt);
                r.RtDeviationMinutes = deviation;
                r.RtDeviationSquared = deviation * deviation;
            }
        }

        // ════════════════════════════════════════════════════════════════
        //  Diagnostic Helpers
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Prints the calibration iteration log as a formatted table to the console.
        /// Used by the benchmark runner (Prompt 5).
        /// </summary>
        public static void PrintCalibrationLog(List<CalibrationIterationLog> log)
        {
            if (log == null || log.Count == 0)
            {
                System.Diagnostics.Debug.WriteLine("  [No calibration iterations recorded]");
                return;
            }

            System.Diagnostics.Debug.WriteLine("  ┌──────────┬──────────┬─────────┬──────────┬──────────┬──────────┬──────────────────┬──────────────────┐");
            System.Diagnostics.Debug.WriteLine("  │ Iter     │ Anchors  │  Slope  │  σ (min) │    R²    │ ± Window │ Model            │ Time             │");
            System.Diagnostics.Debug.WriteLine("  ├──────────┼──────────┼─────────┼──────────┼──────────┼──────────┼──────────────────┼──────────────────┤");

            foreach (var entry in log)
            {
                string modelStr = (entry.ModelType ?? "Linear").PadRight(16);
                string timeStr = $"{entry.TotalTime.TotalMilliseconds:F0} ms".PadRight(16);

                System.Diagnostics.Debug.WriteLine(
                    $"  │ {entry.Iteration,8} │ {entry.AnchorCount,8:N0} │ {entry.Slope,7:F4} │ {entry.SigmaMinutes,8:F3} │ {entry.RSquared,8:F4} │ {entry.WindowHalfWidthMinutes,7:F2}m │ {modelStr} │ {timeStr} │");
            }

            System.Diagnostics.Debug.WriteLine("  └──────────┴──────────┴─────────┴──────────┴──────────┴──────────┴──────────────────┴──────────────────┘");
        }

        /// <summary>
        /// Prints a summary of the pipeline result to the console.
        /// </summary>
        public static void PrintPipelineSummary(PipelineResult result)
        {
            System.Diagnostics.Debug.WriteLine("  ═══ Calibration Pipeline Summary ═══");

            if (result.Calibration != null)
            {
                System.Diagnostics.Debug.WriteLine($"  Final model:     slope={result.Calibration.Slope:F4}  intercept={result.Calibration.Intercept:F2}");
                System.Diagnostics.Debug.WriteLine($"  Final σ:         {result.Calibration.SigmaMinutes:F3} min");
                System.Diagnostics.Debug.WriteLine($"  Final R²:        {result.Calibration.RSquared:F4}");
                System.Diagnostics.Debug.WriteLine($"  Anchor count:    {result.Calibration.AnchorCount:N0}");
                System.Diagnostics.Debug.WriteLine($"  Model reliable:  {result.Calibration.IsReliable}");

                if (result.DetailedCalibration != null &&
                    result.DetailedCalibration.ModelType != RtCalibrationModelType.Linear)
                {
                    System.Diagnostics.Debug.WriteLine($"  Detailed model:  {result.DetailedCalibration.ModelType}");
                }
            }
            else
            {
                System.Diagnostics.Debug.WriteLine("  Calibration FAILED — no model produced.");
            }

            System.Diagnostics.Debug.WriteLine($"  Results:         {result.Results?.Count ?? 0:N0} precursors scored");
            System.Diagnostics.Debug.WriteLine($"  Iterations:      {result.CalibrationLog?.Count ?? 0}");
            System.Diagnostics.Debug.WriteLine($"  Timing:");
            System.Diagnostics.Debug.WriteLine($"    Calibration:   {result.CalibrationTime.TotalSeconds:F1}s");
            System.Diagnostics.Debug.WriteLine($"    Final extract: {result.FinalExtractionTime.TotalSeconds:F1}s");
            System.Diagnostics.Debug.WriteLine($"    RT recalib:    {result.RtDeviationRecalibrationTime.TotalMilliseconds:F0}ms");
            System.Diagnostics.Debug.WriteLine($"    Total:         {result.TotalTime.TotalSeconds:F1}s");
        }

        /// <summary>
        /// Exports calibration convergence data as a TSV string suitable for external plotting.
        /// Each row represents one calibration iteration.
        /// </summary>
        public static string ExportConvergenceTsv(List<CalibrationIterationLog> log)
        {
            if (log == null || log.Count == 0) return string.Empty;

            var sb = new System.Text.StringBuilder();
            sb.AppendLine("Iteration\tAnchorCount\tSlope\tIntercept\tSigma\tRSquared\tWindowHalfWidth\tModelType\tExtractionTimeMs");

            foreach (var entry in log)
            {
                sb.AppendLine(string.Join("\t",
                    entry.Iteration,
                    entry.AnchorCount,
                    entry.Slope.ToString("F4"),
                    entry.Intercept.ToString("F2"),
                    entry.SigmaMinutes.ToString("F3"),
                    entry.RSquared.ToString("F4"),
                    entry.WindowHalfWidthMinutes.ToString("F2"),
                    entry.ModelType ?? "Linear",
                    entry.ExtractionTime.TotalMilliseconds.ToString("F0")));
            }

            return sb.ToString();
        }
    }
}