// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// Phase 15, Prompt 1 + Prompt 4: Diagnostic anchor analysis + RT deviation recalibration
// Placement: MassSpectrometry/Dia/Calibration/RtCalibrationDiagnostics.cs

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MassSpectrometry.Dia.Calibration;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Diagnostic utilities for RT calibration anchor quality analysis
    /// and post-calibration RT deviation correction.
    ///
    /// Prompt 1 (Task 1): DiagnosticAnchorAnalysis — compares ObservedApexRt
    ///   against DIA-NN ground-truth RTs to measure anchor contamination rates.
    ///
    /// Prompt 4 (Task 3): RecalibrateRtDeviations — recomputes RtDeviationMinutes
    ///   using calibrated predictions instead of raw library RT, improving
    ///   target-decoy separation for the classifier.
    ///
    /// These are diagnostic and pipeline helper methods, not production scoring code.
    /// </summary>
    public static class RtCalibrationDiagnostics
    {
        // ════════════════════════════════════════════════════════════════
        //  Prompt 1, Task 1: Diagnostic Anchor Quality Analysis
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Compares ObservedApexRt from our engine against DIA-NN ground truth RTs.
        /// Produces per-anchor: |ObservedApexRt - DiannRt| residuals.
        ///
        /// Reports at multiple ApexScore thresholds (0.5, 0.6, 0.7, 0.8, 0.9):
        ///   - Total anchors, contamination rate (|error| > 0.5, 1.0, 2.0 min)
        ///   - Mean, median, std, P90, P99 of absolute error
        ///   - Text histogram of residuals in 0.1-min bins
        /// </summary>
        /// <param name="broadResults">Results from the broad-window extraction.</param>
        /// <param name="diannRtLookup">
        /// Ground truth observed RT from DIA-NN, keyed by "Sequence_Charge" (e.g., "PEPTIDEK_2").
        /// </param>
        /// <param name="apexScoreThreshold">Primary threshold for detailed analysis.</param>
        /// <param name="output">TextWriter for the diagnostic report.</param>
        public static void DiagnosticAnchorAnalysis(
            List<DiaSearchResult> broadResults,
            Dictionary<string, double> diannRtLookup,
            double apexScoreThreshold,
            TextWriter output)
        {
            if (broadResults == null || diannRtLookup == null || output == null) return;

            output.WriteLine("═══════════════════════════════════════════════════════════");
            output.WriteLine(" RT Calibration Anchor Quality Analysis");
            output.WriteLine("═══════════════════════════════════════════════════════════");
            output.WriteLine();

            // Run at multiple thresholds
            double[] thresholds = { 0.5, 0.6, 0.7, 0.8, 0.9 };

            foreach (double threshold in thresholds)
            {
                var errors = CollectAnchorErrors(broadResults, diannRtLookup, (float)threshold);

                output.WriteLine($"─── ApexScore threshold = {threshold:F1} ───");
                output.WriteLine($"  Matched anchors: {errors.Count:N0}");

                if (errors.Count == 0)
                {
                    output.WriteLine("  (no anchors matched)");
                    output.WriteLine();
                    continue;
                }

                errors.Sort();
                double mean = errors.Average();
                double median = errors[errors.Count / 2];
                double std = Math.Sqrt(errors.Select(e => (e - mean) * (e - mean)).Average());
                double p90 = errors[(int)(errors.Count * 0.90)];
                double p99 = errors[(int)(errors.Count * 0.99)];

                int contam05 = errors.Count(e => e > 0.5);
                int contam10 = errors.Count(e => e > 1.0);
                int contam20 = errors.Count(e => e > 2.0);

                output.WriteLine($"  Mean |error|:    {mean:F3} min");
                output.WriteLine($"  Median |error|:  {median:F3} min");
                output.WriteLine($"  Std |error|:     {std:F3} min");
                output.WriteLine($"  P90 |error|:     {p90:F3} min");
                output.WriteLine($"  P99 |error|:     {p99:F3} min");
                output.WriteLine($"  Contaminated (>0.5 min): {contam05:N0} ({100.0 * contam05 / errors.Count:F1}%)");
                output.WriteLine($"  Contaminated (>1.0 min): {contam10:N0} ({100.0 * contam10 / errors.Count:F1}%)");
                output.WriteLine($"  Contaminated (>2.0 min): {contam20:N0} ({100.0 * contam20 / errors.Count:F1}%)");
                output.WriteLine();

                // Histogram (only for the primary threshold)
                if (Math.Abs(threshold - apexScoreThreshold) < 0.01)
                {
                    PrintHistogram(errors, output, binWidth: 0.1, maxBins: 60);
                }

                output.WriteLine();
            }
        }

        /// <summary>
        /// Collects absolute errors between ObservedApexRt and DIA-NN ground truth RTs
        /// for anchors meeting the given ApexScore threshold.
        /// </summary>
        private static List<double> CollectAnchorErrors(
            List<DiaSearchResult> results,
            Dictionary<string, double> diannRtLookup,
            float minApexScore)
        {
            var errors = new List<double>();

            foreach (var r in results)
            {
                if (r.IsDecoy) continue;
                if (float.IsNaN(r.ApexScore) || r.ApexScore < minApexScore) continue;
                if (float.IsNaN(r.ObservedApexRt)) continue;

                // Key format must match KoinaMspParser.BuildRtLookupFromDiannTsv: "SEQUENCE/charge"
                string key = $"{r.Sequence}/{r.ChargeState}";
                if (!diannRtLookup.TryGetValue(key, out double diannRt)) continue;

                double absError = Math.Abs(r.ObservedApexRt - diannRt);
                errors.Add(absError);
            }

            return errors;
        }

        /// <summary>
        /// Prints a text histogram of error values in fixed-width bins.
        /// </summary>
        private static void PrintHistogram(List<double> values, TextWriter output,
            double binWidth, int maxBins)
        {
            if (values.Count == 0) return;

            double maxVal = values.Max();
            int numBins = Math.Min((int)Math.Ceiling(maxVal / binWidth) + 1, maxBins);

            var counts = new int[numBins];
            foreach (double v in values)
            {
                int bin = Math.Min((int)(v / binWidth), numBins - 1);
                counts[bin]++;
            }

            int maxCount = counts.Max();
            int barWidth = 50;

            output.WriteLine("  Histogram of |ObservedApexRt - DiannRt| (0.1 min bins):");
            for (int b = 0; b < numBins; b++)
            {
                if (counts[b] == 0 && b > 10) continue; // Skip trailing empty bins

                double binStart = b * binWidth;
                int barLen = maxCount > 0 ? (int)Math.Ceiling((double)counts[b] / maxCount * barWidth) : 0;
                string bar = new string('█', barLen);
                output.WriteLine($"  [{binStart,5:F1}, {binStart + binWidth,5:F1}) {counts[b],6:N0} {bar}");
            }
        }

        // ════════════════════════════════════════════════════════════════
        //  Prompt 4, Task 3: RT Deviation Recalibration
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Recomputes RtDeviationMinutes and RtDeviationSquared using calibrated predictions
        /// instead of raw library RT. This must be called after assembly and before feature
        /// extraction to ensure the classifier sees calibrated RT deviations.
        ///
        /// Before calibration:
        ///   RtDeviationMinutes = |ObservedApexRt - LibraryRetentionTime|
        ///   → Similar for targets and decoys (both large and noisy)
        ///
        /// After calibration:
        ///   predictedRt = calibration.ToMinutes(libraryRt)
        ///   RtDeviationMinutes = ObservedApexRt - predictedRt  (signed)
        ///   RtDeviationSquared = RtDeviationMinutes²
        ///   → Targets cluster near 0; decoys spread broadly → better discrimination
        ///
        /// Signed deviation is used (not absolute) because directionality may carry
        /// discriminative information for the classifier.
        /// </summary>
        /// <param name="results">Search results to recalibrate. Modified in place.</param>
        /// <param name="precursors">
        /// Precursor inputs parallel to result indices. Used to retrieve library RT.
        /// </param>
        /// <param name="calibration">The fitted RT calibration model.</param>
        public static void RecalibrateRtDeviations(
            List<DiaSearchResult> results,
            IList<LibraryPrecursorInput> precursors,
            RtCalibrationModel calibration)
        {
            if (results == null || calibration == null)
                return;

            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];

                if (float.IsNaN(r.ObservedApexRt))
                    continue;

                // Use the result's stored LibraryRetentionTime
                double? libRt = r.LibraryRetentionTime;
                if (!libRt.HasValue || double.IsNaN(libRt.Value))
                    continue;

                double predictedRt = calibration.ToMinutes(libRt.Value);

                // Signed deviation
                float deviation = (float)(r.ObservedApexRt - predictedRt);
                r.RtDeviationMinutes = deviation;
                r.RtDeviationSquared = deviation * deviation;
            }
        }

        /// <summary>
        /// Overload accepting <see cref="IRtCalibrationModel"/> for use with piecewise/LOWESS
        /// models that may produce more accurate local predictions than the global linear
        /// approximation from <see cref="RtCalibrationModel"/>.
        /// </summary>
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

                double predictedRt = calibration.ToMinutes(libRt.Value);

                float deviation = (float)(r.ObservedApexRt - predictedRt);
                r.RtDeviationMinutes = deviation;
                r.RtDeviationSquared = deviation * deviation;
            }
        }

        // ════════════════════════════════════════════════════════════════
        //  Additional Prompt 4 Diagnostics: Post-Calibration Analysis
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Prints a summary comparing RT deviations before and after calibration.
        /// Useful for validating that RecalibrateRtDeviations improved target-decoy separation.
        /// </summary>
        /// <param name="results">Results with recalibrated RT deviations.</param>
        /// <param name="output">TextWriter for the report.</param>
        public static void PrintRtDeviationSummary(List<DiaSearchResult> results, TextWriter output)
        {
            if (results == null || output == null) return;

            var targetDevs = new List<double>();
            var decoyDevs = new List<double>();

            foreach (var r in results)
            {
                if (float.IsNaN(r.RtDeviationMinutes)) continue;

                if (r.IsDecoy)
                    decoyDevs.Add(Math.Abs(r.RtDeviationMinutes));
                else
                    targetDevs.Add(Math.Abs(r.RtDeviationMinutes));
            }

            output.WriteLine("─── RT Deviation Summary (post-calibration) ───");

            if (targetDevs.Count > 0)
            {
                targetDevs.Sort();
                double tMean = targetDevs.Average();
                double tMedian = targetDevs[targetDevs.Count / 2];
                double tP90 = targetDevs[(int)(targetDevs.Count * 0.90)];
                output.WriteLine($"  Targets (n={targetDevs.Count:N0}): mean={tMean:F3} median={tMedian:F3} P90={tP90:F3} min");
            }

            if (decoyDevs.Count > 0)
            {
                decoyDevs.Sort();
                double dMean = decoyDevs.Average();
                double dMedian = decoyDevs[decoyDevs.Count / 2];
                double dP90 = decoyDevs[(int)(decoyDevs.Count * 0.90)];
                output.WriteLine($"  Decoys  (n={decoyDevs.Count:N0}): mean={dMean:F3} median={dMedian:F3} P90={dP90:F3} min");
            }

            if (targetDevs.Count > 0 && decoyDevs.Count > 0)
            {
                double separation = decoyDevs.Average() - targetDevs.Average();
                output.WriteLine($"  Separation (decoy_mean - target_mean): {separation:F3} min");
                output.WriteLine($"  (Larger separation → better target-decoy discrimination)");
            }

            output.WriteLine();
        }

        /// <summary>
        /// Computes and reports local σ statistics across the RT range for non-linear models.
        /// Shows how calibration quality varies along the gradient — useful for diagnosing
        /// whether adaptive windows are needed.
        /// </summary>
        /// <param name="detailedModel">The IRtCalibrationModel with GetLocalSigma().</param>
        /// <param name="rtMin">Minimum library RT to sample.</param>
        /// <param name="rtMax">Maximum library RT to sample.</param>
        /// <param name="numBins">Number of RT bins to sample.</param>
        /// <param name="output">TextWriter for the report.</param>
        public static void PrintLocalSigmaProfile(
            IRtCalibrationModel detailedModel,
            double rtMin, double rtMax,
            int numBins,
            TextWriter output)
        {
            if (detailedModel == null || output == null) return;

            output.WriteLine($"─── Local σ Profile ({detailedModel.ModelType}) ───");
            output.WriteLine($"  Global σ = {detailedModel.SigmaMinutes:F3} min");
            output.WriteLine();

            double binWidth = (rtMax - rtMin) / numBins;
            if (binWidth <= 0) return;

            output.WriteLine("  Library RT    Local σ (min)  Ratio to Global");
            output.WriteLine("  ──────────    ─────────────  ───────────────");

            for (int b = 0; b <= numBins; b++)
            {
                double rt = rtMin + b * binWidth;
                double localSigma = detailedModel.GetLocalSigma(rt);
                double ratio = detailedModel.SigmaMinutes > 0
                    ? localSigma / detailedModel.SigmaMinutes
                    : 1.0;

                string bar = ratio > 1.2 ? " ▲" : ratio < 0.8 ? " ▼" : "";
                output.WriteLine($"  {rt,10:F1}    {localSigma,12:F4}  {ratio,14:F2}{bar}");
            }

            output.WriteLine();
        }
    }
}