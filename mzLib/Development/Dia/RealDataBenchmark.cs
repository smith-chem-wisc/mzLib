// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using MassSpectrometry;
using MassSpectrometry.Dia;
using MzLibUtil;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Development.Dia
{
    public static class RealDataBenchmark
    {
        public static void Run(string mzmlPath, string groundTruthTsvPath,
            float ppmTolerance = 20f, float rtToleranceMinutes = 5.0f, int minFragments = 3)
        {
            Console.WriteLine("=== DIA Real Data Benchmark (vs DIA-NN Ground Truth) ===");
            Console.WriteLine();

            // ─── Step 1: Load mzML ──────────────────────────────────────────
            Console.WriteLine($"Loading: {Path.GetFileName(mzmlPath)}");
            var sw = Stopwatch.StartNew();
            var msDataFile = MsDataFileReader.GetDataFile(mzmlPath);
            msDataFile.LoadAllStaticData();
            var allScans = msDataFile.GetAllScansList().ToArray();
            var loadTime = sw.Elapsed;

            Console.WriteLine($"  Load time:    {loadTime.TotalSeconds:F2}s");
            Console.WriteLine($"  Total scans:  {allScans.Length}");
            Console.WriteLine($"  MS1 scans:    {allScans.Count(s => s.MsnOrder == 1)}");
            Console.WriteLine($"  MS2 scans:    {allScans.Count(s => s.MsnOrder == 2)}");
            Console.WriteLine($"  RT range:     {allScans.First().RetentionTime:F2} – {allScans.Last().RetentionTime:F2} min");
            Console.WriteLine();

            // ─── Step 2: Build SoA index ────────────────────────────────────
            sw.Restart();
            using var index = DiaScanIndexBuilder.Build(allScans);
            var buildTime = sw.Elapsed;

            Console.WriteLine("SoA Index built:");
            Console.WriteLine($"  Build time:   {buildTime.TotalMilliseconds:F0} ms");
            Console.WriteLine($"  Windows:      {index.WindowCount}");
            Console.WriteLine($"  MS2 scans:    {index.ScanCount}");
            Console.WriteLine($"  Total peaks:  {index.TotalPeakCount:N0}");
            Console.WriteLine($"  Build rate:   {index.TotalPeakCount / buildTime.TotalSeconds:N0} peaks/sec");

            float indexRtMin = index.GetGlobalRtMin();
            float indexRtMax = index.GetGlobalRtMax();
            Console.WriteLine($"  Index RT:     {indexRtMin:F2} – {indexRtMax:F2} min");
            Console.WriteLine($"  Window boundaries:");
            foreach (int wid in index.GetWindowIds().OrderBy(x => x))
            {
                var (lo, hi) = index.GetWindowBounds(wid);
                index.TryGetScanRangeForWindow(wid, out _, out int count);
                Console.WriteLine($"    Window {wid,2}: {lo:F1} – {hi:F1} m/z ({count} scans)");
            }
            Console.WriteLine();

            // ─── Step 3: Load ground truth ──────────────────────────────────
            Console.WriteLine($"Loading ground truth: {Path.GetFileName(groundTruthTsvPath)}");
            var allPrecursors = LoadGroundTruthTsv(groundTruthTsvPath);
            Console.WriteLine($"  Total precursors: {allPrecursors.Count:N0}");
            Console.WriteLine($"  RT range:         {allPrecursors.Min(p => p.RetentionTime):F2} – " +
                              $"{allPrecursors.Max(p => p.RetentionTime):F2} min");

            float rtFilterMin = indexRtMin - rtToleranceMinutes;
            float rtFilterMax = indexRtMax + rtToleranceMinutes;
            var precursorsInRange = allPrecursors
                .Where(p => p.RetentionTime >= rtFilterMin && p.RetentionTime <= rtFilterMax)
                .ToList();
            Console.WriteLine($"  In RT range:      {precursorsInRange.Count:N0}");

            var precursorsWithWindow = precursorsInRange
                .Where(p => index.FindWindowForPrecursorMz(p.PrecursorMz) >= 0)
                .ToList();
            Console.WriteLine($"  With valid window: {precursorsWithWindow.Count:N0}");
            Console.WriteLine();

            if (precursorsWithWindow.Count == 0) { Console.WriteLine("ERROR: No precursors in range."); return; }

            var dotScorer = new NormalizedDotProductScorer();
            var saScorer = new SpectralAngleScorer();
            int threads = Math.Min(Environment.ProcessorCount, 16);

            // ─── Step 4: Pass 1 — Broad extraction ──────────────────────────
            Console.WriteLine("=== Pass 1: Broad Extraction (for RT calibration) ===");

            var broadParams = new DiaSearchParameters
            {
                PpmTolerance = ppmTolerance,
                RtToleranceMinutes = rtToleranceMinutes,
                MinFragmentsRequired = minFragments,
                ScoringStrategy = ScoringStrategy.ConsensusApex,
            };

            sw.Restart();
            var broadGenResult = DiaLibraryQueryGenerator.Generate(precursorsWithWindow, index, broadParams);
            var broadGenTime = sw.Elapsed;
            Console.WriteLine($"  Query generation:  {broadGenTime.TotalMilliseconds:F0} ms  " +
                              $"({broadGenResult.Queries.Length:N0} queries)");

            sw.Restart();
            ExtractionResult broadExtraction;
            using (var orchestrator = new DiaExtractionOrchestrator(index))
                broadExtraction = orchestrator.ExtractAll(broadGenResult.Queries, maxDegreeOfParallelism: threads);
            var broadExtractTime = sw.Elapsed;
            Console.WriteLine($"  Extraction:        {broadExtractTime.TotalMilliseconds:F0} ms  " +
                              $"({broadExtraction.TotalDataPoints:N0} data points)");

            sw.Restart();
            var broadResults = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                precursorsWithWindow, broadGenResult, broadExtraction, broadParams);
            var broadScoreTime = sw.Elapsed;

            var broadApexScores = broadResults.Where(r => !float.IsNaN(r.ApexDotProductScore))
                .Select(r => r.ApexDotProductScore).ToArray();
            Console.WriteLine($"  Scoring (apex):    {broadScoreTime.TotalMilliseconds:F0} ms  ({broadResults.Count:N0} results)");
            if (broadApexScores.Length > 0)
                Console.WriteLine($"  Apex scores:       median={Median(broadApexScores):F4}, mean={broadApexScores.Average():F4}");
            Console.WriteLine();

            // ─── Step 5: RT Calibration (tighter anchor threshold) ──────────
            Console.WriteLine("=== RT Calibration ===");

            // Use a higher anchor threshold to get cleaner calibration anchors
            float anchorDpThreshold = 0.7f;
            var anchors = new List<(double LibraryRt, double ObservedRt)>();

            for (int g = 0; g < broadGenResult.PrecursorGroups.Length; g++)
            {
                var group = broadGenResult.PrecursorGroups[g];
                var input = precursorsWithWindow[group.InputIndex];
                var matchingResult = broadResults.FirstOrDefault(r =>
                    r.Sequence == input.Sequence && r.ChargeState == input.ChargeState);

                if (matchingResult == null || float.IsNaN(matchingResult.ApexDotProductScore))
                    continue;
                if (matchingResult.ApexDotProductScore < anchorDpThreshold)
                    continue;
                if (!input.RetentionTime.HasValue)
                    continue;

                // Find apex RT from XIC data
                float bestIntensity = 0f;
                float bestRt = (float)input.RetentionTime.Value;
                for (int f = 0; f < group.QueryCount; f++)
                {
                    int qi = group.QueryOffset + f;
                    var fr = broadExtraction.Results[qi];
                    if (fr.DataPointCount == 0) continue;
                    for (int p = 0; p < fr.DataPointCount; p++)
                    {
                        float intensity = broadExtraction.IntensityBuffer[fr.IntensityBufferOffset + p];
                        if (intensity > bestIntensity)
                        {
                            bestIntensity = intensity;
                            bestRt = broadExtraction.RtBuffer[fr.RtBufferOffset + p];
                        }
                    }
                }
                anchors.Add((input.RetentionTime.Value, bestRt));
            }

            Console.WriteLine($"  Anchor candidates (apex DP > {anchorDpThreshold}): {anchors.Count}");

            RtCalibrationModel calibration = null;
            if (anchors.Count >= RtCalibrationModel.MinReliableAnchors)
            {
                var anchorLibRts = anchors.Select(a => a.LibraryRt).ToArray();
                var anchorObsRts = anchors.Select(a => a.ObservedRt).ToArray();

                sw.Restart();
                calibration = RtCalibrationFitter.Fit(
                    (ReadOnlySpan<double>)anchorLibRts, (ReadOnlySpan<double>)anchorObsRts);
                var fitTime = sw.Elapsed;

                if (calibration != null)
                {
                    Console.WriteLine($"  Fit time:     {fitTime.TotalMilliseconds:F1} ms");
                    Console.WriteLine($"  Model:        RT = {calibration.Slope:F4} * libRT + {calibration.Intercept:F4}");
                    Console.WriteLine($"  σ:            {calibration.SigmaMinutes:F4} min");
                    Console.WriteLine($"  R²:           {calibration.RSquared:F4}");
                    Console.WriteLine($"  Anchors used: {calibration.AnchorCount}");
                    Console.WriteLine($"  Reliable:     {(calibration.IsReliable ? "YES" : "NO")}");
                    Console.WriteLine($"  Window at 3σ: ±{calibration.GetMinutesWindowHalfWidth(3.0):F2} min " +
                                      $"(vs ±{rtToleranceMinutes:F1} min broad → " +
                                      $"{rtToleranceMinutes / calibration.GetMinutesWindowHalfWidth(3.0):F1}× narrower)");
                }
            }
            else
            {
                Console.WriteLine($"  Too few anchors ({anchors.Count}) for calibration.");
            }
            Console.WriteLine();

            if (calibration == null || !calibration.IsReliable)
            {
                Console.WriteLine("Calibration not available — skipping detailed analysis.");
                return;
            }

            // ─── Step 6: Sigma Multiplier Sweep ─────────────────────────────
            Console.WriteLine("=== Sigma Multiplier Sweep ===");
            Console.WriteLine($"  Calibration σ = {calibration.SigmaMinutes:F4} min");
            Console.WriteLine();
            Console.WriteLine($"  {"k",-6} {"±Window",-12} {"ExtractMs",-11} {"DataPts",-14} {"Recovery",-10} {"MedApex",-10} {"MedTemporal",-12} {"MedSummed",-10}");
            Console.WriteLine($"  {"─",-6} {"──────",-12} {"────────",-11} {"──────",-14} {"────────",-10} {"───────",-10} {"──────────",-12} {"────────",-10}");

            List<DiaSearchResult> bestResults = null;
            double bestK = 3.0;
            float bestMedianTemporal = 0f;
            ExtractionResult bestExtraction = null;
            DiaLibraryQueryGenerator.GenerationResult bestGenResult = default;

            foreach (double k in new[] { 1.0, 1.5, 2.0, 2.5, 3.0 })
            {
                var kParams = new DiaSearchParameters
                {
                    PpmTolerance = ppmTolerance,
                    RtToleranceMinutes = rtToleranceMinutes,
                    MinFragmentsRequired = minFragments,
                    CalibratedWindowSigmaMultiplier = k,
                    ScoringStrategy = ScoringStrategy.TemporalCosine,
                };
                double halfWidth = calibration.GetMinutesWindowHalfWidth(k);

                sw.Restart();
                var kGenResult = DiaLibraryQueryGenerator.GenerateCalibrated(
                    precursorsWithWindow, index, kParams, calibration);

                ExtractionResult kExtraction;
                using (var orchestrator = new DiaExtractionOrchestrator(index))
                    kExtraction = orchestrator.ExtractAll(kGenResult.Queries, maxDegreeOfParallelism: threads);
                var kExtractTime = sw.Elapsed;

                var kResults = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    precursorsWithWindow, kGenResult, kExtraction, kParams);

                var apexScores = kResults.Where(r => !float.IsNaN(r.ApexDotProductScore)).Select(r => r.ApexDotProductScore).ToArray();
                var temporalScores = kResults.Where(r => !float.IsNaN(r.TemporalCosineScore)).Select(r => r.TemporalCosineScore).ToArray();

                // Also compute summed for comparison
                var summedResults = DiaLibraryQueryGenerator.AssembleResults(
                    precursorsWithWindow, kGenResult, kExtraction.Results, kParams, dotScorer, saScorer);
                var summedScores = summedResults.Where(r => !float.IsNaN(r.DotProductScore)).Select(r => r.DotProductScore).ToArray();

                float recovery = (float)kResults.Count / precursorsWithWindow.Count;
                float medApex = apexScores.Length > 0 ? Median(apexScores) : float.NaN;
                float medTemporal = temporalScores.Length > 0 ? Median(temporalScores) : float.NaN;
                float medSummed = summedScores.Length > 0 ? Median(summedScores) : float.NaN;

                Console.WriteLine($"  {k,-6:F1} {"±" + halfWidth.ToString("F2") + "m",-12} " +
                                  $"{kExtractTime.TotalMilliseconds,-11:F0} {kExtraction.TotalDataPoints,-14:N0} " +
                                  $"{recovery,-10:P1} {medApex,-10:F4} {medTemporal,-12:F4} {medSummed,-10:F4}");

                if (!float.IsNaN(medTemporal) && medTemporal > bestMedianTemporal && recovery >= 0.99f)
                {
                    bestMedianTemporal = medTemporal;
                    bestK = k;
                    bestResults = kResults;
                    bestExtraction = kExtraction;
                    bestGenResult = kGenResult;
                }
            }
            Console.WriteLine();
            Console.WriteLine($"  Best: k={bestK:F1} (median temporal={bestMedianTemporal:F4})");
            Console.WriteLine();

            if (bestResults == null) { Console.WriteLine("No valid results."); return; }

            double bestHalfWidth = calibration.GetMinutesWindowHalfWidth(bestK);

            // ─── Step 7: Nonlinear Power Sweep ──────────────────────────────
            Console.WriteLine($"=== Nonlinear Power Sweep (k={bestK:F1}, ±{bestHalfWidth:F2} min) ===");
            Console.WriteLine();
            Console.WriteLine($"  {"Power",-8} {"Median",-10} {"Mean",-10} {"Q25",-10} {"Q75",-10} {">0.7",-14} {">0.8",-14}");
            Console.WriteLine($"  {"─────",-8} {"──────",-10} {"────",-10} {"───",-10} {"───",-10} {"────",-14} {"────",-14}");

            foreach (float power in new[] { 1.0f, 1.5f, 2.0f, 3.0f })
            {
                var wParams = new DiaSearchParameters
                {
                    PpmTolerance = ppmTolerance,
                    RtToleranceMinutes = rtToleranceMinutes,
                    MinFragmentsRequired = minFragments,
                    CalibratedWindowSigmaMultiplier = bestK,
                    ScoringStrategy = ScoringStrategy.WeightedTemporalCosineWithTransform,
                    NonlinearPower = power,
                };
                var wResults = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    precursorsWithWindow, bestGenResult, bestExtraction, wParams);
                var wScores = wResults.Where(r => !float.IsNaN(r.DotProductScore)).Select(r => r.DotProductScore).ToArray();

                if (wScores.Length > 0)
                {
                    string dp07 = $"{wScores.Count(s => s > 0.7f):N0} ({100.0 * wScores.Count(s => s > 0.7f) / wScores.Length:F1}%)";
                    string dp08 = $"{wScores.Count(s => s > 0.8f):N0} ({100.0 * wScores.Count(s => s > 0.8f) / wScores.Length:F1}%)";
                    Console.WriteLine($"  {power,-8:F1} {Median(wScores),-10:F4} {wScores.Average(),-10:F4} " +
                                      $"{Percentile(wScores, 25),-10:F4} {Percentile(wScores, 75),-10:F4} " +
                                      $"{dp07,-14} {dp08,-14}");
                }
            }
            Console.WriteLine();

            // ─── Step 8: Full Strategy Comparison ────────────────────────────
            Console.WriteLine($"=== Strategy Comparison (k={bestK:F1}, ±{bestHalfWidth:F2} min) ===");
            Console.WriteLine();

            var summedBest = DiaLibraryQueryGenerator.AssembleResults(
                precursorsWithWindow, bestGenResult, bestExtraction.Results,
                new DiaSearchParameters { PpmTolerance = ppmTolerance, MinFragmentsRequired = minFragments },
                dotScorer, saScorer);

            var strategies = new[]
            {
                ("Summed", summedBest.Where(r => !float.IsNaN(r.DotProductScore)).Select(r => r.DotProductScore).ToArray()),
                ("Apex", bestResults.Where(r => !float.IsNaN(r.ApexDotProductScore)).Select(r => r.ApexDotProductScore).ToArray()),
                ("Temporal", bestResults.Where(r => !float.IsNaN(r.TemporalCosineScore)).Select(r => r.TemporalCosineScore).ToArray()),
            };

            Console.WriteLine($"  {"Strategy",-16} {"Median",10} {"Mean",10} {"Q25",10} {"Q75",10} {">0.7",14} {">0.8",14}");
            Console.WriteLine($"  {"--------",-16} {"──────",10} {"────",10} {"───",10} {"───",10} {"────",14} {"────",14}");

            foreach (var (name, scores) in strategies)
            {
                if (scores.Length == 0) continue;
                string dp07 = $"{scores.Count(s => s > 0.7f):N0} ({100.0 * scores.Count(s => s > 0.7f) / scores.Length:F1}%)";
                string dp08 = $"{scores.Count(s => s > 0.8f):N0} ({100.0 * scores.Count(s => s > 0.8f) / scores.Length:F1}%)";
                Console.WriteLine($"  {name,-16} {Median(scores),10:F4} {scores.Average(),10:F4} " +
                                  $"{Percentile(scores, 25),10:F4} {Percentile(scores, 75),10:F4} " +
                                  $"{dp07,14} {dp08,14}");
            }
            Console.WriteLine();

            // ─── Step 9: Hybrid Feature Analysis ────────────────────────────
            Console.WriteLine("=== Hybrid Feature Analysis ===");
            Console.WriteLine("  Every precursor has BOTH apex and temporal scores as independent features.");
            Console.WriteLine();

            var hybrid = bestResults
                .Where(r => !float.IsNaN(r.ApexDotProductScore) && !float.IsNaN(r.TemporalCosineScore))
                .ToList();

            if (hybrid.Count > 0)
            {
                int bothGood = hybrid.Count(r => r.ApexDotProductScore > 0.7f && r.TemporalCosineScore > 0.7f);
                int apexOnly = hybrid.Count(r => r.ApexDotProductScore > 0.7f && r.TemporalCosineScore <= 0.7f);
                int temporalOnly = hybrid.Count(r => r.ApexDotProductScore <= 0.7f && r.TemporalCosineScore > 0.7f);
                int neither = hybrid.Count(r => r.ApexDotProductScore <= 0.7f && r.TemporalCosineScore <= 0.7f);
                int union = bothGood + apexOnly + temporalOnly;

                Console.WriteLine($"  Precursors scored:    {hybrid.Count:N0}");
                Console.WriteLine($"  Both > 0.7:           {bothGood:N0} ({100.0 * bothGood / hybrid.Count:F1}%)");
                Console.WriteLine($"  Apex only > 0.7:      {apexOnly:N0} ({100.0 * apexOnly / hybrid.Count:F1}%)");
                Console.WriteLine($"  Temporal only > 0.7:  {temporalOnly:N0} ({100.0 * temporalOnly / hybrid.Count:F1}%)");
                Console.WriteLine($"  Neither > 0.7:        {neither:N0} ({100.0 * neither / hybrid.Count:F1}%)");
                Console.WriteLine($"  UNION > 0.7:          {union:N0} ({100.0 * union / hybrid.Count:F1}%)");
                Console.WriteLine();

                // Simple combined score: max(apex, temporal)
                var maxScores = hybrid.Select(r => Math.Max(r.ApexDotProductScore, r.TemporalCosineScore)).ToArray();
                var avgScores = hybrid.Select(r => (r.ApexDotProductScore + r.TemporalCosineScore) / 2f).ToArray();
                Console.WriteLine($"  Combined max(apex,temporal): median={Median(maxScores):F4}, " +
                                  $">0.7: {maxScores.Count(s => s > 0.7f):N0} ({100.0 * maxScores.Count(s => s > 0.7f) / maxScores.Length:F1}%)");
                Console.WriteLine($"  Combined avg(apex,temporal): median={Median(avgScores):F4}, " +
                                  $">0.7: {avgScores.Count(s => s > 0.7f):N0} ({100.0 * avgScores.Count(s => s > 0.7f) / avgScores.Length:F1}%)");
                Console.WriteLine();

                Console.WriteLine("  Top 10 by max(apex, temporal):");
                foreach (var r in hybrid.OrderByDescending(r => Math.Max(r.ApexDotProductScore, r.TemporalCosineScore)).Take(10))
                    Console.WriteLine($"    {r.Sequence,-30} z={r.ChargeState} apex={r.ApexDotProductScore:F3} temporal={r.TemporalCosineScore:F3} tpts={r.TimePointsUsed}");
                Console.WriteLine();

                Console.WriteLine("  Largest apex vs temporal disagreements:");
                foreach (var r in hybrid.OrderByDescending(r => MathF.Abs(r.ApexDotProductScore - r.TemporalCosineScore)).Take(10))
                    Console.WriteLine($"    {r.Sequence,-30} z={r.ChargeState} apex={r.ApexDotProductScore:F3} temporal={r.TemporalCosineScore:F3} Δ={MathF.Abs(r.ApexDotProductScore - r.TemporalCosineScore):F3}");
                Console.WriteLine();
            }

            // ─── Step 10: Detailed Correctness ──────────────────────────────
            Console.WriteLine($"=== Correctness Metrics (k={bestK:F1}, ±{bestHalfWidth:F2} min) ===");
            EvaluateResults(bestResults, precursorsWithWindow, bestGenResult);

            // ─── Timing ─────────────────────────────────────────────────────
            Console.WriteLine("=== Timing Summary ===");
            Console.WriteLine($"  mzML load:       {loadTime.TotalSeconds:F2}s");
            Console.WriteLine($"  Index build:     {buildTime.TotalMilliseconds:F0} ms");
            Console.WriteLine($"  Pass 1:          {(broadGenTime + broadExtractTime + broadScoreTime).TotalMilliseconds:F0} ms");
            Console.WriteLine($"  Calibration:     <1 ms");
            Console.WriteLine($"  (sigma sweep and detailed analysis times shown inline above)");
        }

        // ─── Evaluation ─────────────────────────────────────────────────────

        private static void EvaluateResults(
            List<DiaSearchResult> results, List<LibraryPrecursorInput> groundTruth,
            DiaLibraryQueryGenerator.GenerationResult genResult)
        {
            Console.WriteLine();
            Console.WriteLine($"  Precursor Recovery: {results.Count:N0} / {groundTruth.Count:N0} ({100.0 * results.Count / groundTruth.Count:F1}%)");

            if (results.Count == 0) return;

            var detRates = results.Select(r => r.FragmentDetectionRate).ToArray();
            Console.WriteLine($"  Fragment Detection: mean={detRates.Average():P1}, median={Median(detRates):P1}");

            var apex = results.Where(r => !float.IsNaN(r.ApexDotProductScore)).Select(r => r.ApexDotProductScore).ToArray();
            var temporal = results.Where(r => !float.IsNaN(r.TemporalCosineScore)).Select(r => r.TemporalCosineScore).ToArray();

            if (apex.Length > 0)
                Console.WriteLine($"  Apex DP:     median={Median(apex):F4}, Q25={Percentile(apex, 25):F4}, Q75={Percentile(apex, 75):F4}");
            if (temporal.Length > 0)
                Console.WriteLine($"  Temporal DP: median={Median(temporal):F4}, Q25={Percentile(temporal, 25):F4}, Q75={Percentile(temporal, 75):F4}");

            Console.WriteLine();
            Console.WriteLine("  Top 10 by apex score:");
            foreach (var r in results.Where(r => !float.IsNaN(r.ApexDotProductScore)).OrderByDescending(r => r.ApexDotProductScore).Take(10))
                Console.WriteLine($"    {r.Sequence,-30} z={r.ChargeState} apex={r.ApexDotProductScore:F3} temporal={r.TemporalCosineScore:F3} frags={r.FragmentsDetected}/{r.FragmentsQueried}");
            Console.WriteLine();
        }

        // ─── Ground Truth Loader ────────────────────────────────────────────

        public static List<LibraryPrecursorInput> LoadGroundTruthTsv(string path)
        {
            var precursors = new List<LibraryPrecursorInput>();
            bool headerRead = false;
            int colSequence = -1, colPrecursorMz = -1, colCharge = -1, colRt = -1;
            int colIsDecoy = -1, colFragMzs = -1, colFragInts = -1;

            foreach (var line in File.ReadLines(path))
            {
                if (!headerRead)
                {
                    var headers = line.Split('\t');
                    for (int i = 0; i < headers.Length; i++)
                    {
                        switch (headers[i].Trim())
                        {
                            case "Sequence": colSequence = i; break;
                            case "PrecursorMz": colPrecursorMz = i; break;
                            case "Charge": colCharge = i; break;
                            case "RT": colRt = i; break;
                            case "IsDecoy": colIsDecoy = i; break;
                            case "FragmentMzs": colFragMzs = i; break;
                            case "FragmentIntensities": colFragInts = i; break;
                        }
                    }
                    headerRead = true;
                    continue;
                }

                var cols = line.Split('\t');
                if (cols.Length <= colFragInts) continue;

                string sequence = cols[colSequence];
                double precursorMz = double.Parse(cols[colPrecursorMz], CultureInfo.InvariantCulture);
                int charge = int.Parse(cols[colCharge]);
                double rt = double.Parse(cols[colRt], CultureInfo.InvariantCulture);
                bool isDecoy = colIsDecoy >= 0 && cols[colIsDecoy] == "1";

                float[] fragMzs = cols[colFragMzs].Split(';', StringSplitOptions.RemoveEmptyEntries)
                    .Select(s => float.Parse(s, CultureInfo.InvariantCulture)).ToArray();
                float[] fragInts = cols[colFragInts].Split(';', StringSplitOptions.RemoveEmptyEntries)
                    .Select(s => float.Parse(s, CultureInfo.InvariantCulture)).ToArray();

                if (fragMzs.Length == 0 || fragMzs.Length != fragInts.Length) continue;

                precursors.Add(new LibraryPrecursorInput(
                    sequence, precursorMz, charge, rt, isDecoy, fragMzs, fragInts));
            }
            return precursors;
        }

        // ─── Stats ──────────────────────────────────────────────────────────

        private static float Median(float[] values)
        {
            if (values.Length == 0) return float.NaN;
            var sorted = values.OrderBy(x => x).ToArray();
            int mid = sorted.Length / 2;
            return sorted.Length % 2 == 0 ? (sorted[mid - 1] + sorted[mid]) / 2f : sorted[mid];
        }

        private static float Percentile(float[] values, int percentile)
        {
            if (values.Length == 0) return float.NaN;
            var sorted = values.OrderBy(x => x).ToArray();
            int idx = Math.Clamp((int)(sorted.Length * percentile / 100.0), 0, sorted.Length - 1);
            return sorted[idx];
        }
    }
}