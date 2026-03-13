// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
// Location: Development/Dia/Phase12PeakGroupBenchmark.cs

using MassSpectrometry;
using MassSpectrometry.Dia;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace Development.Dia
{
    /// <summary>
    /// Phase 12: Peak Group Detection Benchmark
    /// 
    /// Runs the full DIA pipeline with the Phase 12 peak group detection integrated
    /// into the scoring path. Compares peak-restricted features against full-window
    /// features and measures the impact on FDR yield.
    /// 
    /// Updated Phase 13 Prompt 8: feature vector now has 26 features.
    /// PeakApexScore, MeanFragCorr, MinFragCorr dropped from DiaFeatureVector
    /// but still available on DiaSearchResult for diagnostics.
    /// </summary>
    public static class Phase12PeakGroupBenchmark
    {
        public static void RunAll(
            string rawFilePath,
            string targetMspPath,
            string decoyMspPath,
            string groundTruthTsvPath,
            string outputTsvPath)
        {
            Console.WriteLine("╔══════════════════════════════════════════════════════════════╗");
            Console.WriteLine("║    Phase 12: Peak Group Detection Benchmark                 ║");
            Console.WriteLine("╚══════════════════════════════════════════════════════════════╝");
            Console.WriteLine();

            var totalSw = Stopwatch.StartNew();

            // ════════════════════════════════════════════════════════════
            //  Step 1: Load RT lookup
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 1: Loading RT lookup from ground truth ----------------");
            var sw = Stopwatch.StartNew();

            Dictionary<string, double> rtLookup = null;
            if (!string.IsNullOrEmpty(groundTruthTsvPath) && File.Exists(groundTruthTsvPath))
            {
                rtLookup = KoinaMspParser.BuildRtLookupFromDiannTsv(groundTruthTsvPath);
                Console.WriteLine($"  RT lookup entries: {rtLookup.Count:N0}");
            }
            else
            {
                Console.WriteLine("  WARNING: No ground truth TSV — precursors will lack RT");
            }
            Console.WriteLine($"  Time: {sw.ElapsedMilliseconds}ms");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 2: Load target + decoy libraries
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 2: Loading target + decoy libraries -------------------");
            sw.Restart();

            var targets = KoinaMspParser.Parse(targetMspPath, rtLookup, minIntensity: 0.05f);
            Console.WriteLine($"  Targets: {targets.Count:N0} precursors ({sw.ElapsedMilliseconds}ms)");

            sw.Restart();
            var decoysRaw = KoinaMspParser.Parse(decoyMspPath, rtLookup, minIntensity: 0.05f);
            var decoys = new List<LibraryPrecursorInput>(decoysRaw.Count);
            for (int i = 0; i < decoysRaw.Count; i++)
            {
                var d = decoysRaw[i];
                decoys.Add(new LibraryPrecursorInput(
                    sequence: d.Sequence,
                    precursorMz: d.PrecursorMz,
                    chargeState: d.ChargeState,
                    retentionTime: d.RetentionTime,
                    isDecoy: true,
                    fragmentMzs: d.FragmentMzs,
                    fragmentIntensities: d.FragmentIntensities,
                    irtValue: d.IrtValue));
            }
            Console.WriteLine($"  Decoys:  {decoys.Count:N0} precursors ({sw.ElapsedMilliseconds}ms)");

            var combined = new List<LibraryPrecursorInput>(targets.Count + decoys.Count);
            combined.AddRange(targets);
            combined.AddRange(decoys);
            Console.WriteLine($"  Combined: {combined.Count:N0} precursors");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 3: Load raw file and build index
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 3: Loading raw file and building index ----------------");
            sw.Restart();

            MsDataFile msDataFile;
            string ext = Path.GetExtension(rawFilePath).ToLowerInvariant();
            if (ext == ".raw")
                msDataFile = new ThermoRawFileReader(rawFilePath);
            else
                msDataFile = new Mzml(rawFilePath);
            msDataFile.LoadAllStaticData();
            var loadTime = sw.Elapsed;

            sw.Restart();
            var scans = msDataFile.GetAllScansList().ToArray();
            using var index = DiaScanIndexBuilder.Build(scans);
            var buildTime = sw.Elapsed;

            Console.WriteLine($"  File load: {loadTime.TotalSeconds:F1}s | Index build: {buildTime.TotalMilliseconds:F0}ms");
            Console.WriteLine($"  Scans: {index.ScanCount:N0} | Windows: {index.WindowCount} | Peaks: {index.TotalPeakCount:N0}");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 4: Extraction pipeline
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 4: Running extraction pipeline ------------------------");
            var parameters = new DiaSearchParameters
            {
                PpmTolerance = 20f,
                RtToleranceMinutes = 1.14f,
                MinFragmentsRequired = 3,
                MinScoreThreshold = 0f,
                MaxThreads = -1,
                ScoringStrategy = ScoringStrategy.TemporalCosine,
            };

            sw.Restart();
            var genResult = DiaLibraryQueryGenerator.Generate(combined, index, parameters);
            Console.WriteLine($"  Query generation: {sw.ElapsedMilliseconds}ms | {genResult.Queries.Length:N0} queries");

            sw.Restart();
            using var orchestrator = new DiaExtractionOrchestrator(index);
            var extractionResult = orchestrator.ExtractAll(genResult.Queries,
                maxDegreeOfParallelism: parameters.EffectiveMaxThreads);
            Console.WriteLine($"  Extraction: {sw.ElapsedMilliseconds}ms | {extractionResult.TotalDataPoints:N0} data points");

            // ════════════════════════════════════════════════════════════
            //  Step 5: Assembly + Peak Group Detection (THE KEY CHANGE)
            // ════════════════════════════════════════════════════════════
            Console.WriteLine();
            Console.WriteLine("--- Step 5: Assembly with Peak Group Detection -----------------");
            sw.Restart();
            var results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                combined, genResult, extractionResult, parameters, index);
            var assemblyTime = sw.Elapsed;
            Console.WriteLine($"  Assembly + scoring + peak detection: {assemblyTime.TotalMilliseconds:F0}ms | {results.Count:N0} results");

            // ── Peak group diagnostics ──────────────────────────────────
            PrintPeakGroupDiagnostics(results);

            // ════════════════════════════════════════════════════════════
            //  Step 6: Compute feature vectors
            // ════════════════════════════════════════════════════════════
            Console.WriteLine();
            Console.WriteLine("--- Step 6: Computing feature vectors --------------------------");
            sw.Restart();
            var features = new DiaFeatureVector[results.Count];
            for (int i = 0; i < results.Count; i++)
                features[i] = DiaFeatureExtractor.ComputeFeatures(results[i], i);
            Console.WriteLine($"  Feature computation: {sw.ElapsedMilliseconds}ms | {features.Length:N0} vectors ({DiaFeatureVector.ClassifierFeatureCount} features)");

            // ── Feature distribution comparison ─────────────────────────
            PrintFeatureComparison(results, features);
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 7: Iterative FDR estimation — LDA, NN, GBT
            //  IMPORTANT: RunIterativeFdr writes q-values back to results[]
            //  in-place. Capture counts immediately after each run before
            //  the next classifier overwrites them.
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 7: Iterative FDR estimation (LDA + NN + GBT) ----------");
            Console.WriteLine();

            float[] thresholds = { 0.001f, 0.005f, 0.01f, 0.05f, 0.10f };
            int[] phase12Baseline = { 5548, 9560, 11448, 19514, 27403 };

            // --- LDA ---
            sw.Restart();
            var ldaResult = DiaFdrEngine.RunIterativeFdr(
                results, features,
                DiaClassifierType.LinearDiscriminant);
            long msLda = sw.ElapsedMilliseconds;
            int[] ldaCounts = CountAtThresholds(results, thresholds); // capture before next run
            Console.WriteLine($"  LDA: {msLda}ms | Iterations: {ldaResult.IterationsCompleted} | 1% FDR IDs: {ldaResult.IdentificationsAt1PctFdr:N0}");
            Console.WriteLine();
            Console.WriteLine("--- LDA Per-Iteration Diagnostics ------------------------------");
            Console.WriteLine();
            foreach (var diag in ldaResult.Diagnostics) { DiaFdrEngine.PrintDiagnostics(diag); Console.WriteLine(); }

            // --- NN ---
            sw.Restart();
            var nnResult = DiaFdrEngine.RunIterativeFdr(
                results, features,
                DiaClassifierType.NeuralNetwork);
            long msNn = sw.ElapsedMilliseconds;
            int[] nnCounts = CountAtThresholds(results, thresholds); // capture before next run
            Console.WriteLine($"  NN:  {msNn}ms | Iterations: {nnResult.IterationsCompleted} | 1% FDR IDs: {nnResult.IdentificationsAt1PctFdr:N0}");
            Console.WriteLine();
            Console.WriteLine("--- NN Per-Iteration Diagnostics -------------------------------");
            Console.WriteLine();
            foreach (var diag in nnResult.Diagnostics) { DiaFdrEngine.PrintDiagnostics(diag); Console.WriteLine(); }

            // --- GBT ---
            sw.Restart();
            var gbtResult = DiaFdrEngine.RunIterativeFdr(
                results, features,
                DiaClassifierType.GradientBoostedTree);
            long msGbt = sw.ElapsedMilliseconds;
            int[] gbtCounts = CountAtThresholds(results, thresholds); // GBT is last, no next run
            Console.WriteLine($"  GBT: {msGbt}ms | Iterations: {gbtResult.IterationsCompleted} | 1% FDR IDs: {gbtResult.IdentificationsAt1PctFdr:N0}");
            Console.WriteLine();
            Console.WriteLine("--- GBT Per-Iteration Diagnostics ------------------------------");
            Console.WriteLine();
            foreach (var diag in gbtResult.Diagnostics) { DiaFdrEngine.PrintDiagnostics(diag); Console.WriteLine(); }

            // Best result used for TSV export (use NN q-values — re-run NN to restore)
            var fdrResult = nnResult.IdentificationsAt1PctFdr >= ldaResult.IdentificationsAt1PctFdr
                ? nnResult : ldaResult;

            // ════════════════════════════════════════════════════════════
            //  Step 8: Summary — pre-captured counts per classifier
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- FDR Summary ------------------------------------------------");
            Console.WriteLine();

            Console.WriteLine($"  {"Classifier",-12} | {"q≤0.001",8} | {"q≤0.005",8} | {"q≤0.010",8} | {"q≤0.050",8} | {"q≤0.100",8}");
            Console.WriteLine($"  {"────────────",-12}-+-{"────────",8}-+-{"────────",8}-+-{"────────",8}-+-{"────────",8}-+-{"────────",8}");
            Console.WriteLine($"  {"Ph12 base",-12} | {phase12Baseline[0],8:N0} | {phase12Baseline[1],8:N0} | {phase12Baseline[2],8:N0} | {phase12Baseline[3],8:N0} | {phase12Baseline[4],8:N0}");
            PrintCapturedRow("LDA", ldaCounts, thresholds, phase12Baseline);
            PrintCapturedRow("NN", nnCounts, thresholds, phase12Baseline);
            PrintCapturedRow("GBT", gbtCounts, thresholds, phase12Baseline);
            Console.WriteLine();

            int[] counts = gbtCounts; // counts variable used by TSV export path below

            // ════════════════════════════════════════════════════════════
            //  Step 9: Export TSV
            // ════════════════════════════════════════════════════════════
            if (!string.IsNullOrEmpty(outputTsvPath))
            {
                Console.WriteLine("--- Step 9: Exporting results TSV ------------------------------");
                ExportResultsTsv(results, features, outputTsvPath);
                Console.WriteLine($"  Exported: {outputTsvPath}");
                Console.WriteLine();
            }

            totalSw.Stop();
            Console.WriteLine($"═══ Phase 12 complete. Total time: {totalSw.Elapsed.TotalSeconds:F1}s ═══");
            Console.WriteLine();
        }

        // ════════════════════════════════════════════════════════════════
        //  Peak Group Diagnostics
        // ════════════════════════════════════════════════════════════════

        private static void PrintCapturedRow(
            string name,
            int[] counts,
            float[] thresholds,
            int[] baseline)
        {
            string row = $"  {name,-12}";
            for (int t = 0; t < thresholds.Length; t++)
                row += $" | {counts[t],8:N0}";
            int delta1pct = counts[2] - baseline[2];
            string sign = delta1pct >= 0 ? "+" : "";
            Console.WriteLine($"{row}  ({sign}{delta1pct:N0} vs Ph12)");
        }

        private static int[] CountAtThresholds(List<DiaSearchResult> results, float[] thresholds)
        {
            int[] counts = new int[thresholds.Length];
            for (int i = 0; i < results.Count; i++)
            {
                if (results[i].IsDecoy) continue;
                var fdr = results[i].FdrInfo;
                if (fdr == null) continue;
                for (int t = 0; t < thresholds.Length; t++)
                    if (fdr.QValue <= thresholds[t]) counts[t]++;
            }
            return counts;
        }

        private static void PrintPeakGroupDiagnostics(List<DiaSearchResult> results)
        {
            int totalResults = results.Count;
            int targetResults = 0, decoyResults = 0;
            int targetWithPeak = 0, decoyWithPeak = 0;
            var targetWidths = new List<float>();
            var decoyWidths = new List<float>();
            var targetCandidates = new List<int>();

            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];
                if (r.IsDecoy)
                {
                    decoyResults++;
                    if (r.DetectedPeakGroup.HasValue && r.DetectedPeakGroup.Value.IsValid)
                    {
                        decoyWithPeak++;
                        decoyWidths.Add(r.DetectedPeakGroup.Value.PeakWidthMinutes);
                    }
                }
                else
                {
                    targetResults++;
                    if (r.DetectedPeakGroup.HasValue && r.DetectedPeakGroup.Value.IsValid)
                    {
                        targetWithPeak++;
                        targetWidths.Add(r.DetectedPeakGroup.Value.PeakWidthMinutes);
                        targetCandidates.Add(r.DetectedPeakGroup.Value.CandidateCount);
                    }
                }
            }

            Console.WriteLine();
            Console.WriteLine("  Peak Group Detection:");
            Console.WriteLine($"    Targets: {targetWithPeak:N0}/{targetResults:N0} ({100.0 * targetWithPeak / Math.Max(targetResults, 1):F1}%) have valid peak");
            Console.WriteLine($"    Decoys:  {decoyWithPeak:N0}/{decoyResults:N0} ({100.0 * decoyWithPeak / Math.Max(decoyResults, 1):F1}%) have valid peak");

            if (targetWidths.Count > 0)
            {
                targetWidths.Sort();
                Console.WriteLine($"    Target peak widths: median={targetWidths[targetWidths.Count / 2]:F3} min  " +
                    $"Q25={targetWidths[targetWidths.Count / 4]:F3}  Q75={targetWidths[3 * targetWidths.Count / 4]:F3}");
            }
            if (decoyWidths.Count > 0)
            {
                decoyWidths.Sort();
                Console.WriteLine($"    Decoy peak widths:  median={decoyWidths[decoyWidths.Count / 2]:F3} min  " +
                    $"Q25={decoyWidths[decoyWidths.Count / 4]:F3}  Q75={decoyWidths[3 * decoyWidths.Count / 4]:F3}");
            }
            if (targetCandidates.Count > 0)
            {
                targetCandidates.Sort();
                Console.WriteLine($"    Target candidate peaks: median={targetCandidates[targetCandidates.Count / 2]}  " +
                    $"max={targetCandidates[targetCandidates.Count - 1]}");
            }
        }

        // ════════════════════════════════════════════════════════════════
        //  Feature Comparison (Peak vs Full-Window)
        //  NOTE: PeakApexScore and MeanFragCorr were dropped from the
        //  DiaFeatureVector in Phase 13 Prompt 8 but remain on DiaSearchResult
        //  for diagnostic purposes. This method reads from DiaSearchResult.
        // ════════════════════════════════════════════════════════════════

        private static void PrintFeatureComparison(List<DiaSearchResult> results, DiaFeatureVector[] features)
        {
            Console.WriteLine();
            Console.WriteLine("  Feature comparison (Peak-restricted vs Full-window):");
            Console.WriteLine("  ──────────────────────────────────────────────────────");

            var targetApexFull = new List<float>();
            var targetApexPeak = new List<float>();
            var targetCorrFull = new List<float>();
            var targetCorrPeak = new List<float>();
            var decoyApexPeak = new List<float>();
            var decoyCorrPeak = new List<float>();

            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];
                if (r.IsDecoy)
                {
                    // PeakApexScore read from DiaSearchResult (not in feature vector)
                    if (!float.IsNaN(r.PeakApexScore)) decoyApexPeak.Add(r.PeakApexScore);
                    if (!float.IsNaN(features[i].PeakMeanFragCorr)) decoyCorrPeak.Add(features[i].PeakMeanFragCorr);
                }
                else
                {
                    if (!float.IsNaN(features[i].ApexScore)) targetApexFull.Add(features[i].ApexScore);
                    // PeakApexScore read from DiaSearchResult (not in feature vector)
                    if (!float.IsNaN(r.PeakApexScore)) targetApexPeak.Add(r.PeakApexScore);
                    // MeanFragCorr read from DiaSearchResult (not in feature vector)
                    if (!float.IsNaN(r.MeanFragCorr)) targetCorrFull.Add(r.MeanFragCorr);
                    if (!float.IsNaN(features[i].PeakMeanFragCorr)) targetCorrPeak.Add(features[i].PeakMeanFragCorr);
                }
            }

            targetApexFull.Sort(); targetApexPeak.Sort();
            targetCorrFull.Sort(); targetCorrPeak.Sort();
            decoyApexPeak.Sort(); decoyCorrPeak.Sort();

            if (targetApexFull.Count > 0)
                Console.WriteLine($"    Target ApexScore     full: median={targetApexFull[targetApexFull.Count / 2]:F4}");
            if (targetApexPeak.Count > 0)
                Console.WriteLine($"    Target PeakApexScore peak: median={targetApexPeak[targetApexPeak.Count / 2]:F4}");
            if (decoyApexPeak.Count > 0)
                Console.WriteLine($"    Decoy  PeakApexScore peak: median={decoyApexPeak[decoyApexPeak.Count / 2]:F4}");

            if (targetCorrFull.Count > 0)
                Console.WriteLine($"    Target MeanFragCorr     full: median={targetCorrFull[targetCorrFull.Count / 2]:F4}");
            if (targetCorrPeak.Count > 0)
                Console.WriteLine($"    Target PeakMeanFragCorr peak: median={targetCorrPeak[targetCorrPeak.Count / 2]:F4}");
            if (decoyCorrPeak.Count > 0)
                Console.WriteLine($"    Decoy  PeakMeanFragCorr peak: median={decoyCorrPeak[decoyCorrPeak.Count / 2]:F4}");

            // ── Prompt 3 new features: T/D separation diagnostics ─────────
            Console.WriteLine();
            Console.WriteLine("  New feature diagnostics (Prompt 3):");
            Console.WriteLine("  ──────────────────────────────────────────────────────────────────");
            Console.WriteLine($"  {"Feature",-20} | {"T-mean",8} | {"D-mean",8} | {"Sep (T-D)",10} | {"NaN%",6}");
            Console.WriteLine($"  {new string('-', 20)}-+-{new string('-', 8)}-+-{new string('-', 8)}-+-{new string('-', 10)}-+-{new string('-', 6)}");

            PrintNewFeatureDiag("CoElutionStd", results, features, fv => fv.CoElutionStd, inverted: true);
            PrintNewFeatureDiag("CandidateScoreGap", results, features, fv => fv.CandidateScoreGap, inverted: false);
            Console.WriteLine("  (CoElutionStd: lower=better for targets; CandidateScoreGap: higher=better for targets)");
        }

        private static void PrintNewFeatureDiag(
            string name,
            List<DiaSearchResult> results,
            DiaFeatureVector[] features,
            Func<DiaFeatureVector, float> selector,
            bool inverted)
        {
            double tSum = 0, dSum = 0;
            int tCount = 0, dCount = 0, nanCount = 0;

            for (int i = 0; i < results.Count; i++)
            {
                float v = selector(features[i]);
                if (float.IsNaN(v)) { nanCount++; continue; }
                if (results[i].IsDecoy) { dSum += v; dCount++; }
                else { tSum += v; tCount++; }
            }

            float tMean = tCount > 0 ? (float)(tSum / tCount) : float.NaN;
            float dMean = dCount > 0 ? (float)(dSum / dCount) : float.NaN;
            float sep = (!float.IsNaN(tMean) && !float.IsNaN(dMean)) ? (tMean - dMean) : float.NaN;
            float nanPct = results.Count > 0 ? 100f * nanCount / results.Count : 0f;

            string sepStr = float.IsNaN(sep) ? "   N/A" : $"{sep,+10:F4}";
            string note = !float.IsNaN(sep) && ((inverted && sep < -0.1f) || (!inverted && sep > 0.05f))
                            ? " ✓" : " ⚠";
            Console.WriteLine($"  {name,-20} | {tMean,8:F4} | {dMean,8:F4} | {sepStr} | {nanPct,5:F1}%{note}");
        }

        // ════════════════════════════════════════════════════════════════
        //  TSV Export
        //  Updated Phase 13 Prompt 8: reads dropped features from
        //  DiaSearchResult directly, not from DiaFeatureVector.
        // ════════════════════════════════════════════════════════════════

        private static void ExportResultsTsv(
            List<DiaSearchResult> results,
            DiaFeatureVector[] features,
            string path)
        {
            using var w = new StreamWriter(path);

            w.Write("Sequence\tCharge\tPrecursorMz\tWindowId\tIsDecoy");
            w.Write("\tClassifierScore\tQValue");
            w.Write("\tApexScore\tTemporalScore\tSpectralAngle");
            w.Write("\tMeanFragCorr\tMinFragCorr\tFragDetRate");
            w.Write("\tFragDet\tFragQueried\tTimePointsUsed");
            w.Write("\tObservedApexRt\tLibraryRt\tRtDeviationMin");
            w.Write("\tPeakDetected\tPeakWidth\tPeakSymmetry\tPeakCandidates");
            w.Write("\tPeakApexScore\tPeakMeanFragCorr\tPeakMinFragCorr");
            for (int j = 0; j < DiaFeatureVector.ClassifierFeatureCount; j++)
                w.Write("\tFV_" + DiaFeatureVector.FeatureNames[j]);
            w.WriteLine();

            Span<float> buf = stackalloc float[DiaFeatureVector.ClassifierFeatureCount];
            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];
                double qValue = r.FdrInfo?.QValue ?? 2.0;

                w.Write(r.Sequence);
                w.Write('\t'); w.Write(r.ChargeState);
                w.Write('\t'); w.Write(r.PrecursorMz.ToString("F4"));
                w.Write('\t'); w.Write(r.WindowId);
                w.Write('\t'); w.Write(r.IsDecoy);
                w.Write('\t'); w.Write(r.ClassifierScore.ToString("F6"));
                w.Write('\t'); w.Write(qValue.ToString("F6"));
                w.Write('\t'); w.Write(r.ApexScore.ToString("F4"));
                w.Write('\t'); w.Write(r.TemporalScore.ToString("F4"));
                w.Write('\t'); w.Write(r.SpectralAngleScore.ToString("F4"));
                // MeanFragCorr/MinFragCorr from DiaSearchResult (not in feature vector)
                w.Write('\t'); w.Write(r.MeanFragCorr.ToString("F4"));
                w.Write('\t'); w.Write(r.MinFragCorr.ToString("F4"));
                w.Write('\t'); w.Write(r.FragmentDetectionRate.ToString("F4"));
                w.Write('\t'); w.Write(r.FragmentsDetected);
                w.Write('\t'); w.Write(r.FragmentsQueried);
                w.Write('\t'); w.Write(r.TimePointsUsed);
                w.Write('\t'); w.Write(r.ObservedApexRt.ToString("F4"));
                w.Write('\t'); w.Write(r.LibraryRetentionTime.HasValue
                    ? r.LibraryRetentionTime.Value.ToString("F4") : "NA");
                w.Write('\t'); w.Write(features[i].RtDeviationMinutes.ToString("F4"));

                // Peak group columns — read from DiaSearchResult
                bool hasPeak = r.DetectedPeakGroup.HasValue && r.DetectedPeakGroup.Value.IsValid;
                w.Write('\t'); w.Write(hasPeak);
                w.Write('\t'); w.Write(hasPeak ? r.DetectedPeakGroup.Value.PeakWidthMinutes.ToString("F4") : "NA");
                w.Write('\t'); w.Write(hasPeak ? r.DetectedPeakGroup.Value.SymmetryRatio.ToString("F4") : "NA");
                w.Write('\t'); w.Write(hasPeak ? r.DetectedPeakGroup.Value.CandidateCount.ToString() : "NA");
                // PeakApexScore from DiaSearchResult (dropped from feature vector)
                w.Write('\t'); w.Write(r.PeakApexScore.ToString("F4"));
                w.Write('\t'); w.Write(r.PeakMeanFragCorr.ToString("F4"));
                w.Write('\t'); w.Write(r.PeakMinFragCorr.ToString("F4"));

                features[i].WriteTo(buf);
                for (int j = 0; j < DiaFeatureVector.ClassifierFeatureCount; j++)
                {
                    w.Write('\t');
                    w.Write(buf[j].ToString("G6"));
                }
                w.WriteLine();
            }
        }
    }
}