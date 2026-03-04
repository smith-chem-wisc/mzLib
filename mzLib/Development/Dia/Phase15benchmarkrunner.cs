// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
// Location: Development/Dia/Phase15BenchmarkRunner.cs

using MassSpectrometry;
using MassSpectrometry.Dia;
using MassSpectrometry.Dia.Calibration;
using Readers;
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;  // Used sparingly: .ToArray() in scan loading, library RT anchor collection (not hot paths)

namespace Development.Dia
{
    /// <summary>
    /// Phase 15: RT Calibration Deep Fix — Benchmark Runner and Validation
    ///
    /// Validates the iterative RT calibration system (Prompts 1–4) end-to-end
    /// and produces comparison numbers against the Phase 14 hardcoded baseline.
    ///
    /// Steps:
    ///   1.  Load RT lookup from DIA-NN ground truth TSV
    ///   2.  Load target + decoy libraries from MSP files
    ///   3.  Load raw file and build DiaScanIndex
    ///   4.  Run automatic iterative calibration via DiaCalibrationPipeline
    ///   5.  Print calibration iteration log table
    ///   6.  Diagnostic anchor analysis (ObservedApexRt vs DIA-NN ground truth)
    ///   7.  Feature computation (29 features + recalibrated RT deviation)
    ///   8.  Dead-feature diagnostic table
    ///   9.  LDA FDR — iterative LDA, report 1% FDR IDs
    ///   10. Comparison table: Phase 14 baseline vs Phase 15 automatic calibration
    ///   11. Full FDR threshold table at q ≤ {0.001, 0.005, 0.01, 0.05, 0.10}
    ///   12. Per-step timing summary
    ///   13. TSV export
    ///
    /// A/B Comparison Configurations:
    ///   A: Baseline — Hardcoded ±1.14 min, no iterative calibration (Phase 14 path)
    ///   B: Linear iterative — Iterative calibration, linear model only, starting from ±5.0 min
    ///   C: Non-linear iterative — Iterative calibration, auto model selection, starting from ±5.0 min
    ///
    /// Compilation checklist:
    ///   ✓ Using statements: System, System.Buffers, System.Collections.Generic, System.Diagnostics,
    ///     System.Globalization, System.IO, System.Linq, MassSpectrometry, MassSpectrometry.Dia,
    ///     MassSpectrometry.Dia.Calibration, Readers
    ///   ✓ API consistency: All method calls verified against uploaded files
    ///   ✓ No LINQ in hot paths: Feature loops, FDR counting, TSV export all use explicit for loops
    ///   ✓ Span usage: stackalloc float[29] in diagnostics and TSV export
    ///   ✓ IDisposable: DiaScanIndex (using var), DiaExtractionOrchestrator (using blocks)
    ///   ✓ Pre-run sanity checks: File existence, feature count assertion
    /// </summary>
    public static class Phase15BenchmarkRunner
    {
        // ── Phase 14 baselines for comparison ──────────────────────────
        // These are the q-value counts from Phase 14 LDA at thresholds {0.001, 0.005, 0.01, 0.05, 0.10}
        private static readonly int[] Phase14LdaBaseline = { 8481, 14092, 19873, 24249, 29693 };

        public static void RunAll(
            string rawFilePath,
            string targetMspPath,
            string decoyMspPath,
            string groundTruthTsvPath,
            string outputTsvPath,
            string convergenceTsvPath = null,
            bool runParameterSweep = false)
        {
            Console.WriteLine("╔══════════════════════════════════════════════════════════════╗");
            Console.WriteLine("║  Phase 15: RT Calibration Deep Fix — Benchmark & Validation ║");
            Console.WriteLine("╚══════════════════════════════════════════════════════════════╝");
            Console.WriteLine();

            // ── Pre-run sanity checks ───────────────────────────────────
            if (!File.Exists(rawFilePath))
                throw new FileNotFoundException($"Raw file not found: {rawFilePath}");
            if (!File.Exists(targetMspPath))
                throw new FileNotFoundException($"Target MSP not found: {targetMspPath}");
            if (!File.Exists(decoyMspPath))
                throw new FileNotFoundException($"Decoy MSP not found: {decoyMspPath}");
            if (!string.IsNullOrEmpty(groundTruthTsvPath) && !File.Exists(groundTruthTsvPath))
                Console.WriteLine($"  WARNING: Ground truth TSV not found: {groundTruthTsvPath}");
            if (!string.IsNullOrEmpty(outputTsvPath))
            {
                string outputDir = Path.GetDirectoryName(outputTsvPath);
                if (!string.IsNullOrEmpty(outputDir) && !Directory.Exists(outputDir))
                    Directory.CreateDirectory(outputDir);
            }
            Debug.Assert(DiaFeatureVector.ClassifierFeatureCount == 29,
                $"Expected 29 features, got {DiaFeatureVector.ClassifierFeatureCount}");

            var totalSw = Stopwatch.StartNew();
            var sw = new Stopwatch();

            // Per-step timing
            long msStep1 = 0, msStep2 = 0, msStep3Load = 0, msStep3Index = 0;
            long msStep7 = 0, msStep8 = 0, msStep9 = 0, msStep10 = 0;

            // ════════════════════════════════════════════════════════════
            //  Step 1: Load RT lookup
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 1: Loading RT lookup from ground truth ----------------");
            sw.Restart();

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
            sw.Stop();
            msStep1 = sw.ElapsedMilliseconds;
            Console.WriteLine($"  Time: {msStep1}ms");
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
            sw.Stop();
            msStep2 = sw.ElapsedMilliseconds;
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
            sw.Stop();
            msStep3Load = sw.ElapsedMilliseconds;

            sw.Restart();
            var scans = msDataFile.GetAllScansList().ToArray();
            using var index = DiaScanIndexBuilder.Build(scans);
            sw.Stop();
            msStep3Index = sw.ElapsedMilliseconds;

            Console.WriteLine($"  File load: {msStep3Load / 1000.0:F1}s | Index build: {msStep3Index}ms");
            Console.WriteLine($"  Scans: {index.ScanCount:N0} | Windows: {index.WindowCount} | Peaks: {index.TotalPeakCount:N0}");
            Console.WriteLine();

            int threads = Math.Min(Environment.ProcessorCount, 16);

            // ════════════════════════════════════════════════════════════
            //  Step 4: Run Automatic Iterative Calibration
            //  (Replaces Phase 14 Steps 4-6: broad → calibrate → narrow)
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 4: Automatic Iterative RT Calibration -----------------");
            Console.WriteLine("  Starting from ±5.0 min with auto model selection...");
            Console.WriteLine();
            sw.Restart();

            var calibParams = new DiaSearchParameters
            {
                PpmTolerance = 20f,
                RtToleranceMinutes = 5.0f,    // Starting wide — the whole point of Phase 15
                MinFragmentsRequired = 3,
                MinScoreThreshold = 0f,
                MaxThreads = -1,
                ScoringStrategy = ScoringStrategy.TemporalCosine,
                CalibratedWindowSigmaMultiplier = 3.0,
            };

            // Config C: Non-linear iterative (auto model selection enabled)
            var calibratorC = new IterativeRtCalibrator
            {
                MaxIterations = 4,
                ConvergenceThreshold = 0.05,
                SigmaMultiplier = 3.0,
                InitialTopK = 500,
                InitialApexScoreThreshold = 0.8f,
                RefinedApexScoreThreshold = 0.5f,
                MinWindowHalfWidthMinutes = 0.3,
                EnableNonLinearModelSelection = true,
            };

            DiaCalibrationPipeline.PipelineResult pipelineResult;
            using (var orchestrator = new DiaExtractionOrchestrator(index))
            {
                pipelineResult = DiaCalibrationPipeline.RunWithAutomaticCalibration(
                    combined, index, calibParams, orchestrator, calibratorC);
            }

            sw.Stop();
            var results = pipelineResult.Results;
            var calibration = pipelineResult.Calibration;
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 5: Print Calibration Iteration Log
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 5: Calibration Iteration Log --------------------------");
            Console.WriteLine();

            DiaCalibrationPipeline.PrintCalibrationLog(pipelineResult.CalibrationLog);
            Console.WriteLine();
            DiaCalibrationPipeline.PrintPipelineSummary(pipelineResult);
            Console.WriteLine();
            // ── Soft Assertions (Regression Tests) ──────────────────────
            if (calibration != null)
            {
                if (calibration.SigmaMinutes >= 1.0)
                    Console.WriteLine($"  WARN: Calibration σ too large: {calibration.SigmaMinutes:F3} min (expected < 1.0)");
                else
                    Console.WriteLine($"  PASS: Calibration σ = {calibration.SigmaMinutes:F3} min (< 1.0)");
                if (Math.Abs(calibration.Slope - 1.0) >= 0.1)
                    Console.WriteLine($"  WARN: Calibration slope far from 1.0: {calibration.Slope:F4} (expected ~0.99 for same-instrument library)");
                else
                    Console.WriteLine($"  PASS: Calibration slope = {calibration.Slope:F4} (within 0.1 of 1.0)");
                if (pipelineResult.CalibrationLog != null && pipelineResult.CalibrationLog.Count > 5)
                    Console.WriteLine($"  WARN: Calibration required {pipelineResult.CalibrationLog.Count} iterations (expected ≤ 5)");
                else
                    Console.WriteLine($"  PASS: Calibration converged in {pipelineResult.CalibrationLog?.Count ?? 0} iterations (≤ 5)");
            }
            else
            {
                Console.WriteLine("  FAIL: Calibration produced no model!");
            }
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 6: Diagnostic Anchor Analysis
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 6: Diagnostic Anchor Analysis -------------------------");
            Console.WriteLine();

            if (rtLookup != null && rtLookup.Count > 0 && results != null)
            {
                RtCalibrationDiagnostics.DiagnosticAnchorAnalysis(
                    results, rtLookup, apexScoreThreshold: 0.5, output: Console.Out);

                // Also print RT deviation summary (post-calibration)
                RtCalibrationDiagnostics.PrintRtDeviationSummary(results, Console.Out);

                // If we have a detailed model, print local σ profile
                if (pipelineResult.DetailedCalibration != null)
                {
                    float globalRtMin = index.GetGlobalRtMin();
                    float globalRtMax = index.GetGlobalRtMax();
                    RtCalibrationDiagnostics.PrintLocalSigmaProfile(
                        pipelineResult.DetailedCalibration,
                        globalRtMin, globalRtMax, numBins: 10,
                        output: Console.Out);
                }
            }
            else
            {
                Console.WriteLine("  Skipped (no ground truth RT lookup available).");
            }
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 7: Feature computation (29 features)
            //  Note: RT deviations already recalibrated by the pipeline
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 7: Computing feature vectors --------------------------");
            sw.Restart();

            if (DiaFeatureVector.ClassifierFeatureCount != 29)
                Console.WriteLine($"  WARNING: Expected 29 features, got {DiaFeatureVector.ClassifierFeatureCount}");

            var features = new DiaFeatureVector[results.Count];
            for (int i = 0; i < results.Count; i++)
                features[i] = DiaFeatureExtractor.ComputeFeatures(results[i], i);
            sw.Stop();
            msStep7 = sw.ElapsedMilliseconds;
            Console.WriteLine($"  Feature computation: {msStep7}ms | {features.Length:N0} vectors ({DiaFeatureVector.ClassifierFeatureCount} features)");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 8: Dead-feature diagnostic table
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 8: Dead Feature Fix Verification ----------------------");
            sw.Restart();

            PrintDeadFeatureDiagnostics(results, features);

            sw.Stop();
            msStep8 = sw.ElapsedMilliseconds;
            Console.WriteLine($"  Time: {msStep8}ms");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 9: LDA FDR
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 9: LDA FDR --------------------------------------------");
            sw.Restart();

            var ldaResult = DiaFdrEngine.RunIterativeFdr(
                results, features,
                classifierType: DiaClassifierType.LinearDiscriminant);

            sw.Stop();
            msStep9 = sw.ElapsedMilliseconds;
            Console.WriteLine($"  FDR estimation: {msStep9}ms");
            Console.WriteLine($"  Iterations: {ldaResult.IterationsCompleted}");
            Console.WriteLine($"  Final 1% FDR IDs: {ldaResult.IdentificationsAt1PctFdr:N0}");
            Console.WriteLine();

            // Capture q-value counts at multiple thresholds
            float[] thresholds = { 0.001f, 0.005f, 0.01f, 0.05f, 0.10f };
            int[] configCCounts = CountIdsAtThresholds(results, thresholds);
            int configCAt1Pct = configCCounts[2]; // q <= 0.01

            // Print per-iteration diagnostics
            Console.WriteLine("  Per-Iteration Diagnostics:");
            foreach (var diag in ldaResult.Diagnostics)
            {
                DiaFdrEngine.PrintDiagnostics(diag);
                Console.WriteLine();
            }

            // ── Soft assertion: ID count ────────────────────────────────
            Console.WriteLine(configCAt1Pct >= 19000
                ? $"  PASS: {configCAt1Pct:N0} IDs ≥ 19,000 target"
                : $"  WARN: {configCAt1Pct:N0} IDs < 19,000 target (Phase 14 baseline: 19,873)");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 10: A/B/C Comparison Table
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 10: A/B/C Configuration Comparison --------------------");
            Console.WriteLine();

            // Config A: Baseline (Phase 14 hardcoded ±1.14 min)
            Console.WriteLine("  Running Config A: Baseline (±1.14 min, no iterative calibration)...");
            var configAResult = RunBaselineConfig(combined, index, threads);
            int[] configACounts = null;
            string configASummary = "N/A";

            if (configAResult.Results != null && configAResult.Results.Count > 0)
            {
                var featA = new DiaFeatureVector[configAResult.Results.Count];
                for (int i = 0; i < configAResult.Results.Count; i++)
                    featA[i] = DiaFeatureExtractor.ComputeFeatures(configAResult.Results[i], i);

                DiaFdrEngine.RunIterativeFdr(
                    configAResult.Results, featA,
                    classifierType: DiaClassifierType.LinearDiscriminant);
                configACounts = CountIdsAtThresholds(configAResult.Results, thresholds);
                configASummary = $"slope=N/A σ=N/A window=±1.14m";
            }
            Console.WriteLine($"    Config A done: {(configACounts != null ? configACounts[2].ToString("N0") : "?")} IDs at 1% FDR");

            // Config B: Linear-only iterative
            Console.WriteLine("  Running Config B: Linear iterative (±5.0 min start, linear only)...");
            var configBResult = RunLinearOnlyConfig(combined, index, calibParams, threads);
            int[] configBCounts = null;
            string configBSummary = "N/A";

            if (configBResult.Results != null && configBResult.Results.Count > 0)
            {
                var featB = new DiaFeatureVector[configBResult.Results.Count];
                for (int i = 0; i < configBResult.Results.Count; i++)
                    featB[i] = DiaFeatureExtractor.ComputeFeatures(configBResult.Results[i], i);

                DiaFdrEngine.RunIterativeFdr(
                    configBResult.Results, featB,
                    classifierType: DiaClassifierType.LinearDiscriminant);
                configBCounts = CountIdsAtThresholds(configBResult.Results, thresholds);

                if (configBResult.Calibration != null)
                {
                    configBSummary = $"slope={configBResult.Calibration.Slope:F4} " +
                                     $"σ={configBResult.Calibration.SigmaMinutes:F3} " +
                                     $"window=±{configBResult.Calibration.GetMinutesWindowHalfWidth(3.0):F2}m";
                }
            }
            Console.WriteLine($"    Config B done: {(configBCounts != null ? configBCounts[2].ToString("N0") : "?")} IDs at 1% FDR");

            // Config C is the primary run (already computed above)
            string configCSummary = "N/A";
            if (calibration != null)
            {
                string modelType = pipelineResult.DetailedCalibration?.ModelType.ToString() ?? "Linear";
                configCSummary = $"slope={calibration.Slope:F4} " +
                                 $"σ={calibration.SigmaMinutes:F3} " +
                                 $"window=±{calibration.GetMinutesWindowHalfWidth(3.0):F2}m " +
                                 $"model={modelType}";
            }

            Console.WriteLine();
            Console.WriteLine("  ┌────────┬────────────────────────────────────────────────────────────────────────────┐");
            Console.WriteLine("  │ Config │ Description                                                              │");
            Console.WriteLine("  ├────────┼────────────────────────────────────────────────────────────────────────────┤");
            Console.WriteLine("  │   A    │ Baseline: hardcoded ±1.14 min, no iterative calibration (Phase 14 path)  │");
            Console.WriteLine("  │   B    │ Linear iterative: ±5.0 min start, linear model only                     │");
            Console.WriteLine("  │   C    │ Non-linear iterative: ±5.0 min start, auto model selection              │");
            Console.WriteLine("  └────────┴────────────────────────────────────────────────────────────────────────────┘");
            Console.WriteLine();

            Console.WriteLine("  Config  |  q≤0.001  q≤0.005  q≤0.01   q≤0.05   q≤0.10   | Calibration");
            Console.WriteLine("  ─────────────────────────────────────────────────────────────────────────────");

            PrintConfigRow("A: Baseline", configACounts, Phase14LdaBaseline, configASummary);
            PrintConfigRow("B: Linear", configBCounts, Phase14LdaBaseline, configBSummary);
            PrintConfigRow("C: AutoModel", configCCounts, Phase14LdaBaseline, configCSummary);
            PrintConfigRow("P14 Ref", Phase14LdaBaseline, null, "Phase 14 LDA baseline");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 11: Full FDR Threshold Table (Config C)
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 11: Full FDR Table (Config C: Auto Model) -------------");
            Console.WriteLine();

            Console.WriteLine($"  Q-value threshold  |  Target IDs  |  vs Phase 14 LDA");
            Console.WriteLine($"  ─────────────────────────────────────────────────────");
            for (int t = 0; t < thresholds.Length; t++)
            {
                int delta = configCCounts[t] - Phase14LdaBaseline[t];
                string sign = delta >= 0 ? "+" : "";
                Console.WriteLine($"  q ≤ {thresholds[t]:F3}           |  {configCCounts[t],6:N0}       |  {sign}{delta:N0}");
            }
            Console.WriteLine();

            // Classifier score statistics
            double targetScoreSum = 0, decoyScoreSum = 0;
            int targetCount = 0, decoyCount = 0;
            for (int i = 0; i < results.Count; i++)
            {
                float cs = results[i].ClassifierScore;
                if (float.IsNaN(cs)) continue;
                if (results[i].IsDecoy) { decoyScoreSum += cs; decoyCount++; }
                else { targetScoreSum += cs; targetCount++; }
            }
            float meanTarget = targetCount > 0 ? (float)(targetScoreSum / targetCount) : 0f;
            float meanDecoy = decoyCount > 0 ? (float)(decoyScoreSum / decoyCount) : 0f;
            Console.WriteLine($"  Total scored: {results.Count:N0} (targets={targetCount:N0}, decoys={decoyCount:N0})");
            Console.WriteLine($"  Mean classifier score: target={meanTarget:F4}  decoy={meanDecoy:F4}  gap={meanTarget - meanDecoy:F4}");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 12: Per-step timing summary
            // ════════════════════════════════════════════════════════════
            totalSw.Stop();
            Console.WriteLine("--- Step 12: Timing Summary ------------------------------------");
            Console.WriteLine();
            Console.WriteLine($"  Step 1   RT lookup load           {msStep1,8}ms");
            Console.WriteLine($"  Step 2   Library load             {msStep2,8}ms");
            Console.WriteLine($"  Step 3   Raw file load            {msStep3Load,8}ms");
            Console.WriteLine($"  Step 3   Index build              {msStep3Index,8}ms");
            Console.WriteLine($"  Step 4   Calibration pipeline     {(long)pipelineResult.TotalTime.TotalMilliseconds,8}ms");
            Console.WriteLine($"    ├─ Calibration loop             {(long)pipelineResult.CalibrationTime.TotalMilliseconds,8}ms");
            Console.WriteLine($"    ├─ Final extraction             {(long)pipelineResult.FinalExtractionTime.TotalMilliseconds,8}ms");
            Console.WriteLine($"    └─ RT deviation recalibration   {(long)pipelineResult.RtDeviationRecalibrationTime.TotalMilliseconds,8}ms");
            Console.WriteLine($"  Step 7   Feature computation      {msStep7,8}ms");
            Console.WriteLine($"  Step 8   Dead feature diag        {msStep8,8}ms");
            Console.WriteLine($"  Step 9   LDA FDR                  {msStep9,8}ms");
            Console.WriteLine($"  ──────────────────────────────────────────────");
            Console.WriteLine($"  TOTAL                             {totalSw.ElapsedMilliseconds,8}ms ({totalSw.Elapsed.TotalSeconds:F1}s)");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 13: TSV Export + Convergence TSV
            // ════════════════════════════════════════════════════════════
            if (!string.IsNullOrEmpty(outputTsvPath))
            {
                Console.WriteLine("--- Step 13: Exporting results TSV -----------------------------");
                ExportResultsTsv(results, features, outputTsvPath, "LDA");
                Console.WriteLine($"  Exported: {outputTsvPath}");
                Console.WriteLine($"  Rows: {results.Count:N0} | Columns: {5 + 4 + DiaFeatureVector.ClassifierFeatureCount + 6 + 16}");
                Console.WriteLine();
            }

            // Calibration convergence TSV
            if (!string.IsNullOrEmpty(convergenceTsvPath) && pipelineResult.CalibrationLog != null)
            {
                Console.WriteLine("--- Step 13b: Exporting calibration convergence TSV ------------");
                string convergenceTsv = DiaCalibrationPipeline.ExportConvergenceTsv(pipelineResult.CalibrationLog);
                File.WriteAllText(convergenceTsvPath, convergenceTsv);
                Console.WriteLine($"  Exported: {convergenceTsvPath}");
                Console.WriteLine();
            }

            // ════════════════════════════════════════════════════════════
            //  Optional: Parameter Sensitivity Sweep
            // ════════════════════════════════════════════════════════════
            if (runParameterSweep)
            {
                Console.WriteLine("--- Parameter Sensitivity Sweep --------------------------------");
                Console.WriteLine();
                RunParameterSensitivitySweep(combined, index, calibParams, threads);
                Console.WriteLine();
            }

            Console.WriteLine($"═══ Phase 15 complete. Total time: {totalSw.Elapsed.TotalSeconds:F1}s ═══");
            Console.WriteLine();
        }

        // ════════════════════════════════════════════════════════════════
        //  Config A: Baseline (Phase 14 path — hardcoded ±1.14 min)
        // ════════════════════════════════════════════════════════════════

        private readonly struct ConfigResult
        {
            public readonly List<DiaSearchResult> Results;
            public readonly RtCalibrationModel Calibration;
            public readonly List<CalibrationIterationLog> CalibrationLog;

            public ConfigResult(List<DiaSearchResult> results,
                RtCalibrationModel calibration,
                List<CalibrationIterationLog> log)
            {
                Results = results;
                Calibration = calibration;
                CalibrationLog = log;
            }
        }

        private static ConfigResult RunBaselineConfig(
            IList<LibraryPrecursorInput> combined,
            DiaScanIndex index,
            int threads)
        {
            var baselineParams = new DiaSearchParameters
            {
                PpmTolerance = 20f,
                RtToleranceMinutes = 1.14f,    // Phase 14 hardcoded value
                MinFragmentsRequired = 3,
                MinScoreThreshold = 0f,
                MaxThreads = -1,
                ScoringStrategy = ScoringStrategy.TemporalCosine,
            };

            DiaLibraryQueryGenerator.GenerationResult genResult;
            ExtractionResult extraction;
            List<DiaSearchResult> results;

            genResult = DiaLibraryQueryGenerator.Generate(combined, index, baselineParams);

            using (var orchestrator = new DiaExtractionOrchestrator(index))
                extraction = orchestrator.ExtractAll(genResult.Queries,
                    maxDegreeOfParallelism: threads);

            results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                combined, genResult, extraction, baselineParams, index);

            // Calibrate (same as Phase 14 Step 5)
            float anchorDpThreshold = 0.5f;
            var anchors = new List<(double LibraryRt, double ObservedRt)>();
            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];
                if (r.IsDecoy) continue;
                if (float.IsNaN(r.ApexScore) || r.ApexScore < anchorDpThreshold) continue;
                if (!r.LibraryRetentionTime.HasValue) continue;
                if (float.IsNaN(r.ObservedApexRt)) continue;
                anchors.Add((r.LibraryRetentionTime.Value, r.ObservedApexRt));
            }

            RtCalibrationModel calibA = null;
            if (anchors.Count >= RtCalibrationModel.MinReliableAnchors)
            {
                var anchorLibRts = anchors.Select(a => a.LibraryRt).ToArray();
                var anchorObsRts = anchors.Select(a => a.ObservedRt).ToArray();
                calibA = RtCalibrationFitter.Fit(
                    (ReadOnlySpan<double>)anchorLibRts, (ReadOnlySpan<double>)anchorObsRts);
            }

            // If calibration available, do calibrated re-extraction
            if (calibA != null && calibA.IsReliable)
            {
                var calibBaseParams = new DiaSearchParameters
                {
                    PpmTolerance = 20f,
                    RtToleranceMinutes = 5.0f,
                    MinFragmentsRequired = 3,
                    MinScoreThreshold = 0f,
                    MaxThreads = -1,
                    ScoringStrategy = ScoringStrategy.TemporalCosine,
                    CalibratedWindowSigmaMultiplier = 2.0,
                };

                genResult = DiaLibraryQueryGenerator.GenerateCalibrated(
                    combined, index, calibBaseParams, calibA);

                using (var orchestrator = new DiaExtractionOrchestrator(index))
                    extraction = orchestrator.ExtractAll(genResult.Queries,
                        maxDegreeOfParallelism: threads);

                results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    combined, genResult, extraction, calibBaseParams, index);
            }

            return new ConfigResult(results, calibA, null);
        }

        // ════════════════════════════════════════════════════════════════
        //  Config B: Linear-only iterative calibration
        // ════════════════════════════════════════════════════════════════

        private static ConfigResult RunLinearOnlyConfig(
            IList<LibraryPrecursorInput> combined,
            DiaScanIndex index,
            DiaSearchParameters calibParams,
            int threads)
        {
            var calibratorB = new IterativeRtCalibrator
            {
                MaxIterations = 4,
                ConvergenceThreshold = 0.05,
                SigmaMultiplier = 3.0,
                InitialTopK = 500,
                InitialApexScoreThreshold = 0.8f,
                RefinedApexScoreThreshold = 0.5f,
                MinWindowHalfWidthMinutes = 0.3,
                EnableNonLinearModelSelection = false,    // Linear only
            };

            DiaCalibrationPipeline.PipelineResult result;
            using (var orchestrator = new DiaExtractionOrchestrator(index))
            {
                result = DiaCalibrationPipeline.RunWithAutomaticCalibration(
                    combined, index, calibParams, orchestrator, calibratorB);
            }

            return new ConfigResult(result.Results, result.Calibration, result.CalibrationLog);
        }

        // ════════════════════════════════════════════════════════════════
        //  Parameter Sensitivity Sweep (Task 5)
        // ════════════════════════════════════════════════════════════════

        private static void RunParameterSensitivitySweep(
            IList<LibraryPrecursorInput> combined,
            DiaScanIndex index,
            DiaSearchParameters baseParams,
            int threads)
        {
            float[] qThresholds = { 0.001f, 0.005f, 0.01f, 0.05f, 0.10f };

            Console.WriteLine("  ┌─────────────────────────────────┬────────┬─────────┬──────────┐");
            Console.WriteLine("  │ Parameter Variation             │ IDs@1% │  σ(min) │ Iters    │");
            Console.WriteLine("  ├─────────────────────────────────┼────────┼─────────┼──────────┤");

            // Sweep InitialTopK
            int[] topKValues = { 200, 500, 1000, 2000 };
            foreach (int topK in topKValues)
            {
                var result = RunSweepConfig(combined, index, baseParams, threads,
                    cal => cal.InitialTopK = topK);
                PrintSweepRow($"InitialTopK={topK}", result, qThresholds);
            }

            // Sweep InitialApexScoreThreshold
            float[] apexThresholds = { 0.7f, 0.8f, 0.9f };
            foreach (float thresh in apexThresholds)
            {
                var result = RunSweepConfig(combined, index, baseParams, threads,
                    cal => cal.InitialApexScoreThreshold = thresh);
                PrintSweepRow($"InitApexThresh={thresh:F1}", result, qThresholds);
            }

            // Sweep SigmaMultiplier
            double[] sigmaMultipliers = { 2.0, 2.5, 3.0 };
            foreach (double mult in sigmaMultipliers)
            {
                var result = RunSweepConfig(combined, index, baseParams, threads,
                    cal => cal.SigmaMultiplier = mult);
                PrintSweepRow($"SigmaMultiplier={mult:F1}", result, qThresholds);
            }

            // Sweep ConvergenceThreshold
            double[] convergenceThresholds = { 0.02, 0.05, 0.10 };
            foreach (double thresh in convergenceThresholds)
            {
                var result = RunSweepConfig(combined, index, baseParams, threads,
                    cal => cal.ConvergenceThreshold = thresh);
                PrintSweepRow($"ConvThresh={thresh:F2}", result, qThresholds);
            }

            // Sweep MinWindowHalfWidthMinutes
            double[] minWindows = { 0.2, 0.3, 0.5 };
            foreach (double minW in minWindows)
            {
                var result = RunSweepConfig(combined, index, baseParams, threads,
                    cal => cal.MinWindowHalfWidthMinutes = minW);
                PrintSweepRow($"MinWindow={minW:F1}m", result, qThresholds);
            }

            Console.WriteLine("  └─────────────────────────────────┴────────┴─────────┴──────────┘");
        }

        private readonly struct SweepResult
        {
            public readonly int IdsAt1Pct;
            public readonly double Sigma;
            public readonly int Iterations;

            public SweepResult(int ids, double sigma, int iterations)
            {
                IdsAt1Pct = ids;
                Sigma = sigma;
                Iterations = iterations;
            }
        }

        private static SweepResult RunSweepConfig(
            IList<LibraryPrecursorInput> combined,
            DiaScanIndex index,
            DiaSearchParameters baseParams,
            int threads,
            Action<IterativeRtCalibrator> configure)
        {
            var calibrator = new IterativeRtCalibrator
            {
                MaxIterations = 4,
                ConvergenceThreshold = 0.05,
                SigmaMultiplier = 3.0,
                InitialTopK = 500,
                InitialApexScoreThreshold = 0.8f,
                RefinedApexScoreThreshold = 0.5f,
                MinWindowHalfWidthMinutes = 0.3,
                EnableNonLinearModelSelection = true,
            };
            configure(calibrator);

            DiaCalibrationPipeline.PipelineResult result;
            using (var orchestrator = new DiaExtractionOrchestrator(index))
            {
                result = DiaCalibrationPipeline.RunWithAutomaticCalibration(
                    combined, index, baseParams, orchestrator, calibrator);
            }

            if (result.Results == null || result.Results.Count == 0)
                return new SweepResult(0, double.NaN, result.CalibrationLog?.Count ?? 0);

            var feat = new DiaFeatureVector[result.Results.Count];
            for (int i = 0; i < result.Results.Count; i++)
                feat[i] = DiaFeatureExtractor.ComputeFeatures(result.Results[i], i);

            var fdrResult = DiaFdrEngine.RunIterativeFdr(
                result.Results, feat,
                classifierType: DiaClassifierType.LinearDiscriminant);

            double sigma = result.Calibration?.SigmaMinutes ?? double.NaN;
            return new SweepResult(fdrResult.IdentificationsAt1PctFdr, sigma,
                result.CalibrationLog?.Count ?? 0);
        }

        private static void PrintSweepRow(string label, SweepResult result, float[] thresholds)
        {
            Console.WriteLine(
                $"  │ {label,-33} │ {result.IdsAt1Pct,6:N0} │ {result.Sigma,7:F3} │ {result.Iterations,8} │");
        }

        // ════════════════════════════════════════════════════════════════
        //  Dead Feature Diagnostics (reused from Phase 14)
        // ════════════════════════════════════════════════════════════════

        private static void PrintDeadFeatureDiagnostics(
            List<DiaSearchResult> results, DiaFeatureVector[] features)
        {
            var deadFeatures = new[]
            {
                (Name: "PeakWidth",          Index: 5,  Phase13Value: "0.0 all"),
                (Name: "TimePointsUsed",     Index: 12, Phase13Value: "0.0 all"),
                (Name: "RtDeviationMinutes", Index: 10, Phase13Value: "5.0 all"),
                (Name: "RtDeviationSquared", Index: 11, Phase13Value: "25.0 all"),
            };

            Console.WriteLine();
            Console.WriteLine("  Phase 15 Dead Feature Fix Verification:");
            Console.WriteLine("  ─────────────────────────────────────────────────────────────────────────");
            Console.WriteLine("  Feature              | Phase 13   | P15 T-mean  | P15 D-mean  | Sep     | NaN%   | Fixed?");
            Console.WriteLine("  ─────────────────────────────────────────────────────────────────────────");

            Span<float> buf = stackalloc float[DiaFeatureVector.ClassifierFeatureCount];

            foreach (var (name, idx, p13val) in deadFeatures)
            {
                double targetSum = 0, decoySum = 0;
                int targetN = 0, decoyN = 0;
                int nanCount = 0;

                for (int i = 0; i < results.Count; i++)
                {
                    features[i].WriteTo(buf);
                    float val = buf[idx];

                    if (float.IsNaN(val))
                    {
                        nanCount++;
                        continue;
                    }

                    if (results[i].IsDecoy)
                    {
                        decoySum += val;
                        decoyN++;
                    }
                    else
                    {
                        targetSum += val;
                        targetN++;
                    }
                }

                float tMean = targetN > 0 ? (float)(targetSum / targetN) : float.NaN;
                float dMean = decoyN > 0 ? (float)(decoySum / decoyN) : float.NaN;
                float separation = (!float.IsNaN(tMean) && !float.IsNaN(dMean))
                    ? MathF.Abs(tMean - dMean) : 0f;
                float nanPct = 100f * nanCount / results.Count;
                bool fixed_ = separation > 0.001f && nanPct < 50f;

                Console.WriteLine($"  {name,-22} | {p13val,-10} | {tMean,11:F4} | {dMean,11:F4} | {separation,7:F4} | {nanPct,5:F1}% | {(fixed_ ? "YES" : "NO")}");
            }

            Console.WriteLine("  ─────────────────────────────────────────────────────────────────────────");
            Console.WriteLine();

            // Full 29-feature summary
            Console.WriteLine("  Full 29-Feature Summary (target vs decoy means):");
            Console.WriteLine("  ───────────────────────────────────────────────────────────");

            for (int f = 0; f < DiaFeatureVector.ClassifierFeatureCount; f++)
            {
                double tSum = 0, dSum = 0;
                int tN = 0, dN = 0;

                for (int i = 0; i < results.Count; i++)
                {
                    features[i].WriteTo(buf);
                    float val = buf[f];
                    if (float.IsNaN(val)) continue;
                    if (results[i].IsDecoy) { dSum += val; dN++; }
                    else { tSum += val; tN++; }
                }

                float tm = tN > 0 ? (float)(tSum / tN) : float.NaN;
                float dm = dN > 0 ? (float)(dSum / dN) : float.NaN;
                float sep = (!float.IsNaN(tm) && !float.IsNaN(dm)) ? MathF.Abs(tm - dm) : 0f;

                Console.WriteLine($"    [{f,2}] {DiaFeatureVector.FeatureNames[f],-26} T={tm,9:F4}  D={dm,9:F4}  sep={sep,7:F4}");
            }
            Console.WriteLine();
        }

        // ════════════════════════════════════════════════════════════════
        //  Comparison Table Helpers
        // ════════════════════════════════════════════════════════════════

        private static void PrintConfigRow(string name, int[] counts, int[] baseline, string summary)
        {
            string row = $"  {name,-14} |";
            if (counts != null)
            {
                for (int t = 0; t < counts.Length; t++)
                {
                    if (baseline != null)
                    {
                        int delta = counts[t] - baseline[t];
                        string sign = delta >= 0 ? "+" : "";
                        row += $"  {counts[t],5:N0}({sign}{delta:N0})";
                    }
                    else
                    {
                        row += $"  {counts[t],8:N0}";
                    }
                }
            }
            else
            {
                row += "     N/A     N/A     N/A     N/A     N/A";
            }
            row += $"  | {summary}";
            Console.WriteLine(row);
        }

        // ════════════════════════════════════════════════════════════════
        //  Count IDs at q-value thresholds
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

        // ════════════════════════════════════════════════════════════════
        //  TSV Export (same as Phase 14 — 60 columns)
        // ════════════════════════════════════════════════════════════════

        private static void ExportResultsTsv(
            List<DiaSearchResult> results,
            DiaFeatureVector[] features,
            string path,
            string classifierType)
        {
            using var w = new StreamWriter(path);
            var inv = CultureInfo.InvariantCulture;

            // ── Header ──────────────────────────────────────────────────
            // Group 1: Identification
            w.Write("Sequence\tCharge\tPrecursorMz\tWindowId\tIsDecoy");
            // Group 2: Classifier
            w.Write("\tClassifierScore\tQValue\tPeptideQValue\tClassifierType");
            // Group 3: All 29 feature vector values
            for (int j = 0; j < DiaFeatureVector.ClassifierFeatureCount; j++)
                w.Write("\tFV_" + DiaFeatureVector.FeatureNames[j]);
            // Group 4: Raw result properties for dead-feature debugging
            w.Write("\tObservedApexRt\tLibraryRT\tRtDeviationMinutes_Raw");
            w.Write("\tPeakWidth_Raw\tTimePointsUsed_Raw\tHasPeakGroup");
            // Group 5: Phase 13 feature properties from DiaSearchResult
            w.Write("\tMeanMassErrorPpm\tMassErrorStdPpm\tMaxAbsMassErrorPpm");
            w.Write("\tBestFragCorrelationSum\tMedianFragRefCorr\tMinFragRefCorr\tStdFragRefCorr");
            w.Write("\tMeanSignalRatioDev\tMaxSignalRatioDev\tStdSignalRatioDev");
            w.Write("\tSmoothedMeanFragCorr\tSmoothedMinFragCorr\tLog2SignalToNoise");
            w.Write("\tBestFragWeightedCosine\tBoundarySignalRatio\tApexToMeanRatio");
            w.WriteLine();

            // ── Data rows (explicit for loop — no LINQ) ─────────────────
            Span<float> buf = stackalloc float[DiaFeatureVector.ClassifierFeatureCount];
            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];
                double qValue = r.FdrInfo?.QValue ?? 2.0;
                string pepQValue = r.FdrInfo?.PeptideQValue.HasValue == true
                    ? r.FdrInfo.PeptideQValue.Value.ToString("F6", inv) : "NA";

                // Group 1: Identification
                w.Write(r.Sequence);
                w.Write('\t'); w.Write(r.ChargeState.ToString(inv));
                w.Write('\t'); w.Write(r.PrecursorMz.ToString("F4", inv));
                w.Write('\t'); w.Write(r.WindowId.ToString(inv));
                w.Write('\t'); w.Write(r.IsDecoy ? "True" : "False");

                // Group 2: Classifier
                w.Write('\t'); w.Write(FloatToTsv(r.ClassifierScore, inv));
                w.Write('\t'); w.Write(qValue.ToString("F6", inv));
                w.Write('\t'); w.Write(pepQValue);
                w.Write('\t'); w.Write(classifierType);

                // Group 3: Feature vector values
                features[i].WriteTo(buf);
                for (int j = 0; j < DiaFeatureVector.ClassifierFeatureCount; j++)
                {
                    w.Write('\t');
                    w.Write(FloatToTsv(buf[j], inv));
                }

                // Group 4: Raw result properties for dead-feature debugging
                w.Write('\t'); w.Write(FloatToTsv(r.ObservedApexRt, inv));
                w.Write('\t'); w.Write(r.LibraryRetentionTime.HasValue
                    ? r.LibraryRetentionTime.Value.ToString("F4", inv) : "NA");
                w.Write('\t'); w.Write(FloatToTsv(r.RtDeviationMinutes, inv));
                w.Write('\t'); w.Write(FloatToTsv(r.PeakWidth, inv));
                w.Write('\t'); w.Write(r.TimePointsUsed.ToString(inv));
                bool hasPeak = r.DetectedPeakGroup.HasValue && r.DetectedPeakGroup.Value.IsValid;
                w.Write('\t'); w.Write(hasPeak ? "True" : "False");

                // Group 5: Phase 13 feature properties
                w.Write('\t'); w.Write(FloatToTsv(r.MeanMassErrorPpm, inv));
                w.Write('\t'); w.Write(FloatToTsv(r.MassErrorStdPpm, inv));
                w.Write('\t'); w.Write(FloatToTsv(r.MaxAbsMassErrorPpm, inv));
                w.Write('\t'); w.Write(FloatToTsv(r.BestFragCorrelationSum, inv));
                w.Write('\t'); w.Write(FloatToTsv(r.MedianFragRefCorr, inv));
                w.Write('\t'); w.Write(FloatToTsv(r.MinFragRefCorr, inv));
                w.Write('\t'); w.Write(FloatToTsv(r.StdFragRefCorr, inv));
                w.Write('\t'); w.Write(FloatToTsv(r.MeanSignalRatioDev, inv));
                w.Write('\t'); w.Write(FloatToTsv(r.MaxSignalRatioDev, inv));
                w.Write('\t'); w.Write(FloatToTsv(r.StdSignalRatioDev, inv));
                w.Write('\t'); w.Write(FloatToTsv(r.SmoothedMeanFragCorr, inv));
                w.Write('\t'); w.Write(FloatToTsv(r.SmoothedMinFragCorr, inv));
                w.Write('\t'); w.Write(FloatToTsv(r.Log2SignalToNoise, inv));
                w.Write('\t'); w.Write(FloatToTsv(r.BestFragWeightedCosine, inv));
                w.Write('\t'); w.Write(FloatToTsv(r.BoundarySignalRatio, inv));
                w.Write('\t'); w.Write(FloatToTsv(r.ApexToMeanRatio, inv));

                w.WriteLine();
            }
        }

        /// <summary>Formats a float for TSV: "NA" for NaN, InvariantCulture otherwise.</summary>
        private static string FloatToTsv(float value, IFormatProvider fmt) =>
            float.IsNaN(value) ? "NA" : value.ToString("G6", fmt);
    }
}