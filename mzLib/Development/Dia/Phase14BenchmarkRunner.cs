// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
// Location: Development/Dia/Phase14BenchmarkRunner.cs

using MassSpectrometry;
using MassSpectrometry.Dia;
using Readers;
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;  // Used sparingly: .ToArray() in scan loading and RT anchor collection (not hot paths)

namespace Development.Dia
{
    /// <summary>
    /// Phase 14: Dead Feature Fixes + GBT Hyperparameter Tuning Benchmark
    /// 
    /// Steps 1-7:  Standard pipeline (load → index → broad extract → calibrate → calibrated extract → assemble → features)
    /// Step 8:     Dead-feature diagnostic table (verify PeakWidth, TimePointsUsed, RtDeviationMinutes, RtDeviationSquared)
    /// Step 9:     LDA baseline FDR
    /// Step 10:    GBT sweep (configs A–E)
    /// Step 11:    Comparison table (Phase 13 LDA, Phase 14 LDA, GBT configs)
    /// Step 12:    Full FDR threshold table for best classifier
    /// Step 13:    Per-step timing summary
    /// Step 14:    TSV export with full column set (finalized in Prompt 4)
    /// 
    /// Compilation checklist (Prompt 4 Part B):
    ///   ✓ Using statements: System, System.Buffers, System.Collections.Generic, System.Diagnostics,
    ///     System.Globalization, System.IO, System.Linq, MassSpectrometry, MassSpectrometry.Dia, Readers
    ///   ✓ API consistency: All method calls verified against uploaded files
    ///   ✓ No LINQ in hot paths: Feature loops, FDR counting, TSV export all use explicit for loops
    ///   ✓ Span usage: stackalloc float[29] in diagnostics and TSV export
    ///   ✓ IDisposable: DiaScanIndex (using var), DiaExtractionOrchestrator (using blocks)
    ///   ✓ Pre-run sanity checks: File existence, feature count assertion
    /// </summary>
    public static class Phase14BenchmarkRunner
    {
        public static void RunAll(
            string rawFilePath,
            string targetMspPath,
            string decoyMspPath,
            string groundTruthTsvPath,
            string outputTsvPath)
        {
            Console.WriteLine("╔══════════════════════════════════════════════════════════════╗");
            Console.WriteLine("║    Phase 14: Dead Feature Fixes + GBT Tuning Benchmark      ║");
            Console.WriteLine("╚══════════════════════════════════════════════════════════════╝");
            Console.WriteLine();

            // ── Pre-run sanity checks (Prompt 4 Part C) ─────────────────
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
            long msStep4Gen = 0, msStep4Extract = 0, msStep4Assembly = 0;
            long msStep5 = 0;
            long msStep6Gen = 0, msStep6Extract = 0, msStep6Assembly = 0;
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
            //  Step 4: Broad extraction pass (for RT calibration)
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 4: Broad extraction pass (±5.0 min) -------------------");

            var broadParams = new DiaSearchParameters
            {
                PpmTolerance = 20f,
                RtToleranceMinutes = 1.14f,
                MinFragmentsRequired = 3,
                MinScoreThreshold = 0f,
                MaxThreads = -1,
                ScoringStrategy = ScoringStrategy.TemporalCosine,
            };

            sw.Restart();
            var broadGenResult = DiaLibraryQueryGenerator.Generate(combined, index, broadParams);
            sw.Stop();
            msStep4Gen = sw.ElapsedMilliseconds;
            Console.WriteLine($"  Query generation: {msStep4Gen}ms | {broadGenResult.Queries.Length:N0} queries");

            sw.Restart();
            ExtractionResult broadExtraction;
            using (var orchestrator = new DiaExtractionOrchestrator(index))
                broadExtraction = orchestrator.ExtractAll(broadGenResult.Queries,
                    maxDegreeOfParallelism: threads);
            sw.Stop();
            msStep4Extract = sw.ElapsedMilliseconds;
            Console.WriteLine($"  Extraction: {msStep4Extract}ms | {broadExtraction.TotalDataPoints:N0} data points");

            sw.Restart();
            var broadResults = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                combined, broadGenResult, broadExtraction, broadParams, index);
            sw.Stop();
            msStep4Assembly = sw.ElapsedMilliseconds;
            Console.WriteLine($"  Assembly + scoring: {msStep4Assembly}ms | {broadResults.Count:N0} results");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 5: RT Calibration
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 5: RT Calibration -------------------------------------");
            sw.Restart();

            float anchorDpThreshold = 0.5f;
            var anchors = new List<(double LibraryRt, double ObservedRt)>();

            for (int i = 0; i < broadResults.Count; i++)
            {
                var r = broadResults[i];
                if (r.IsDecoy) continue;
                if (float.IsNaN(r.ApexScore) || r.ApexScore < anchorDpThreshold) continue;
                if (!r.LibraryRetentionTime.HasValue) continue;
                if (float.IsNaN(r.ObservedApexRt)) continue;
                anchors.Add((r.LibraryRetentionTime.Value, r.ObservedApexRt));
            }

            Console.WriteLine($"  Anchor candidates (apex ≥ {anchorDpThreshold}): {anchors.Count}");

            RtCalibrationModel calibration = null;
            if (anchors.Count >= RtCalibrationModel.MinReliableAnchors)
            {
                // LINQ usage here is deliberate: one-time array creation, not a hot path
                var anchorLibRts = anchors.Select(a => a.LibraryRt).ToArray();
                var anchorObsRts = anchors.Select(a => a.ObservedRt).ToArray();

                calibration = RtCalibrationFitter.Fit(
                    (ReadOnlySpan<double>)anchorLibRts, (ReadOnlySpan<double>)anchorObsRts);

                if (calibration != null)
                {
                    Console.WriteLine($"  Model: RT = {calibration.Slope:F4} * libRT + {calibration.Intercept:F4}");
                    Console.WriteLine($"  σ = {calibration.SigmaMinutes:F4} min | R² = {calibration.RSquared:F4}");
                    Console.WriteLine($"  Anchors used: {calibration.AnchorCount} | Reliable: {(calibration.IsReliable ? "YES" : "NO")}");
                    Console.WriteLine($"  Window at 2σ: ±{calibration.GetMinutesWindowHalfWidth(2.0):F2} min");
                }
            }
            else
            {
                Console.WriteLine($"  Too few anchors ({anchors.Count}) for calibration.");
            }
            sw.Stop();
            msStep5 = sw.ElapsedMilliseconds;
            Console.WriteLine($"  Time: {msStep5}ms");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 6: Calibrated extraction pass (±2σ window)
            // ════════════════════════════════════════════════════════════
            List<DiaSearchResult> results;
            DiaLibraryQueryGenerator.GenerationResult genResult;
            ExtractionResult extraction;

            if (calibration != null && calibration.IsReliable)
            {
                Console.WriteLine("--- Step 6: Calibrated extraction pass (±2σ) -------------------");

                var calibParams = new DiaSearchParameters
                {
                    PpmTolerance = 20f,
                    RtToleranceMinutes = 5.0f,
                    MinFragmentsRequired = 3,
                    MinScoreThreshold = 0f,
                    MaxThreads = -1,
                    ScoringStrategy = ScoringStrategy.TemporalCosine,
                    CalibratedWindowSigmaMultiplier = 2.0,
                };

                sw.Restart();
                genResult = DiaLibraryQueryGenerator.GenerateCalibrated(
                    combined, index, calibParams, calibration);
                sw.Stop();
                msStep6Gen = sw.ElapsedMilliseconds;
                Console.WriteLine($"  Query generation: {msStep6Gen}ms | {genResult.Queries.Length:N0} queries");

                sw.Restart();
                using (var orchestrator = new DiaExtractionOrchestrator(index))
                    extraction = orchestrator.ExtractAll(genResult.Queries,
                        maxDegreeOfParallelism: threads);
                sw.Stop();
                msStep6Extract = sw.ElapsedMilliseconds;
                Console.WriteLine($"  Extraction: {msStep6Extract}ms | {extraction.TotalDataPoints:N0} data points");

                sw.Restart();
                results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    combined, genResult, extraction, calibParams, index);
                sw.Stop();
                msStep6Assembly = sw.ElapsedMilliseconds;
                Console.WriteLine($"  Assembly + scoring: {msStep6Assembly}ms | {results.Count:N0} results");
            }
            else
            {
                Console.WriteLine("--- Step 6: SKIPPED (calibration unavailable) — using broad pass");
                results = broadResults;
            }
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 7: Feature computation (29 features)
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
            //  Step 9: LDA baseline FDR
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 9: LDA Baseline FDR -----------------------------------");
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

            // Capture LDA q-value counts at multiple thresholds
            float[] thresholds = { 0.001f, 0.005f, 0.01f, 0.05f, 0.10f };
            int[] ldaCounts = CountIdsAtThresholds(results, thresholds);
            int ldaAt1Pct = ldaCounts[2]; // q <= 0.01

            // Print per-iteration diagnostics
            Console.WriteLine("  Per-Iteration Diagnostics:");
            foreach (var diag in ldaResult.Diagnostics)
            {
                DiaFdrEngine.PrintDiagnostics(diag);
                Console.WriteLine();
            }

            // ════════════════════════════════════════════════════════════
            //  Step 10: GBT Sweep
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 10: GBT Hyperparameter Sweep --------------------------");
            sw.Restart();

            var gbtConfigs = Phase14GbtSweep.GetStandardConfigs();
            var sweepResults = Phase14GbtSweep.RunGbtSweep(results, features);

            sw.Stop();
            msStep10 = sw.ElapsedMilliseconds;
            Console.WriteLine($"  GBT sweep: {msStep10}ms ({gbtConfigs.Length} configs)");
            Console.WriteLine();

            // Print sweep summary
            Phase14GbtSweep.PrintSummaryTable(sweepResults, ldaCounts);

            // ════════════════════════════════════════════════════════════
            //  Step 11: Comparison table
            // ════════════════════════════════════════════════════════════
            Console.WriteLine();
            Console.WriteLine("--- Step 11: Phase 13 → Phase 14 Comparison --------------------");
            Console.WriteLine();

            // Phase 13 baselines
            int[] phase13Lda = { 8481, 14092, 16357, 24249, 29693 };

            Console.WriteLine("  Classifier          |  q≤0.001  q≤0.005  q≤0.01   q≤0.05   q≤0.10");
            Console.WriteLine("  ────────────────────────────────────────────────────────────────────");
            PrintComparisonRow("Phase 13 LDA", phase13Lda, null);
            PrintComparisonRow("Phase 14 LDA", ldaCounts, phase13Lda);

            // Find best GBT config
            int bestGbtIdx = -1;
            int bestGbtAt1Pct = 0;
            for (int i = 0; i < sweepResults.Length; i++)
            {
                if (sweepResults[i].IdsAtQ001 > bestGbtAt1Pct)
                {
                    bestGbtAt1Pct = sweepResults[i].IdsAtQ001;
                    bestGbtIdx = i;
                }
            }

            for (int i = 0; i < sweepResults.Length; i++)
            {
                var sr = sweepResults[i];
                int[] gbtCounts = { sr.IdsAtQ0001, sr.IdsAtQ0005, sr.IdsAtQ001, sr.IdsAtQ005, sr.IdsAtQ010 };
                string marker = (i == bestGbtIdx) ? " ★" : "";
                PrintComparisonRow($"GBT {gbtConfigs[i].Name}{marker}", gbtCounts, phase13Lda);
            }
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 12: Full FDR threshold table for best classifier
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 12: Full FDR Table (Best Classifier) ------------------");
            Console.WriteLine();

            // Determine best overall classifier
            string bestClassifierName;
            int bestOverallAt1Pct;

            if (bestGbtIdx >= 0 && bestGbtAt1Pct > ldaAt1Pct)
            {
                bestClassifierName = $"GBT {gbtConfigs[bestGbtIdx].Name}";
                bestOverallAt1Pct = bestGbtAt1Pct;

                // Re-run best GBT config on the actual results so they have the right scores
                Console.WriteLine($"  Best classifier: {bestClassifierName} ({bestOverallAt1Pct:N0} IDs at 1% FDR)");
                Console.WriteLine($"  Re-applying best GBT config to results for export...");
                DiaFdrEngine.RunIterativeFdr(results, features, gbtParams: gbtConfigs[bestGbtIdx]);
            }
            else
            {
                bestClassifierName = "LDA";
                bestOverallAt1Pct = ldaAt1Pct;
                Console.WriteLine($"  Best classifier: LDA ({bestOverallAt1Pct:N0} IDs at 1% FDR)");

                // Re-run LDA so results have LDA scores for export
                DiaFdrEngine.RunIterativeFdr(results, features,
                    classifierType: DiaClassifierType.LinearDiscriminant);
            }

            int[] bestCounts = CountIdsAtThresholds(results, thresholds);
            Console.WriteLine();
            Console.WriteLine($"  Q-value threshold  |  Target IDs  |  vs Phase 13 LDA");
            Console.WriteLine($"  ─────────────────────────────────────────────────────");
            for (int t = 0; t < thresholds.Length; t++)
            {
                int delta = bestCounts[t] - phase13Lda[t];
                string sign = delta >= 0 ? "+" : "";
                Console.WriteLine($"  q ≤ {thresholds[t]:F3}           |  {bestCounts[t],6:N0}       |  {sign}{delta:N0}");
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
            //  Step 13: Per-step timing summary
            // ════════════════════════════════════════════════════════════
            totalSw.Stop();
            Console.WriteLine("--- Step 13: Timing Summary ------------------------------------");
            Console.WriteLine();
            Console.WriteLine($"  Step 1  RT lookup load        {msStep1,8}ms");
            Console.WriteLine($"  Step 2  Library load          {msStep2,8}ms");
            Console.WriteLine($"  Step 3  Raw file load         {msStep3Load,8}ms");
            Console.WriteLine($"  Step 3  Index build           {msStep3Index,8}ms");
            Console.WriteLine($"  Step 4  Broad query gen       {msStep4Gen,8}ms");
            Console.WriteLine($"  Step 4  Broad extraction      {msStep4Extract,8}ms");
            Console.WriteLine($"  Step 4  Broad assembly        {msStep4Assembly,8}ms");
            Console.WriteLine($"  Step 5  RT calibration        {msStep5,8}ms");
            Console.WriteLine($"  Step 6  Calib query gen       {msStep6Gen,8}ms");
            Console.WriteLine($"  Step 6  Calib extraction      {msStep6Extract,8}ms");
            Console.WriteLine($"  Step 6  Calib assembly        {msStep6Assembly,8}ms");
            Console.WriteLine($"  Step 7  Feature computation   {msStep7,8}ms");
            Console.WriteLine($"  Step 8  Dead feature diag     {msStep8,8}ms");
            Console.WriteLine($"  Step 9  LDA FDR               {msStep9,8}ms");
            Console.WriteLine($"  Step 10 GBT sweep             {msStep10,8}ms");
            Console.WriteLine($"  ────────────────────────────────────────────");
            Console.WriteLine($"  TOTAL                         {totalSw.ElapsedMilliseconds,8}ms ({totalSw.Elapsed.TotalSeconds:F1}s)");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 14: TSV export (finalized — Prompt 4)
            //
            //  Column groups:
            //    1. Identification:  Sequence, Charge, PrecursorMz, WindowId, IsDecoy
            //    2. Classifier:      ClassifierScore, QValue, PeptideQValue, ClassifierType
            //    3. Feature vector:  FV_0 .. FV_28  (29 features via WriteTo)
            //    4. Raw properties:  ObservedApexRt, LibraryRT, RtDeviationMinutes_Raw,
            //                        PeakWidth_Raw, TimePointsUsed_Raw, HasPeakGroup
            //    5. Phase 13 props:  16 properties from DiaSearchResult
            // ════════════════════════════════════════════════════════════
            if (!string.IsNullOrEmpty(outputTsvPath))
            {
                Console.WriteLine("--- Step 14: Exporting results TSV -----------------------------");
                ExportResultsTsv(results, features, outputTsvPath, bestClassifierName);
                Console.WriteLine($"  Exported: {outputTsvPath}");
                Console.WriteLine($"  Rows: {results.Count:N0} | Columns: {5 + 4 + DiaFeatureVector.ClassifierFeatureCount + 6 + 16}");
                Console.WriteLine();
            }

            Console.WriteLine($"═══ Phase 14 complete. Total time: {totalSw.Elapsed.TotalSeconds:F1}s ═══");
            Console.WriteLine();
        }

        // ════════════════════════════════════════════════════════════════
        //  Step 8: Dead Feature Diagnostics
        // ════════════════════════════════════════════════════════════════

        private static void PrintDeadFeatureDiagnostics(
            List<DiaSearchResult> results, DiaFeatureVector[] features)
        {
            // Feature indices for the 4 previously-dead features
            // [5] = PeakWidth, [10] = RtDeviationMinutes, [11] = RtDeviationSquared, [12] = TimePointsUsed
            var deadFeatures = new[]
            {
                (Name: "PeakWidth",          Index: 5,  Phase13Value: "0.0 all"),
                (Name: "TimePointsUsed",     Index: 12, Phase13Value: "0.0 all"),
                (Name: "RtDeviationMinutes", Index: 10, Phase13Value: "5.0 all"),
                (Name: "RtDeviationSquared", Index: 11, Phase13Value: "25.0 all"),
            };

            Console.WriteLine();
            Console.WriteLine("  Phase 14 Dead Feature Fix Verification:");
            Console.WriteLine("  ─────────────────────────────────────────────────────────────────────────");
            Console.WriteLine("  Feature              | Phase 13   | P14 T-mean  | P14 D-mean  | Sep     | NaN%   | Fixed?");
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

            // Also print full 29-feature summary for reference
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
        //  Comparison table helper
        // ════════════════════════════════════════════════════════════════

        private static void PrintComparisonRow(string name, int[] counts, int[] baseline)
        {
            string row = $"  {name,-20} |";
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
            Console.WriteLine(row);
        }

        // ════════════════════════════════════════════════════════════════
        //  TSV Export (finalized — Prompt 4 Part A)
        //
        //  Column layout verified against DiaSearchResult, DiaFeatureVector,
        //  DiaFdrInfo, and Prompt 4 specification:
        //
        //  Group 1 — Identification (5 cols):
        //    Sequence, Charge, PrecursorMz, WindowId, IsDecoy
        //
        //  Group 2 — Classifier (4 cols):
        //    ClassifierScore, QValue, PeptideQValue, ClassifierType
        //
        //  Group 3 — Feature vector (29 cols):
        //    FV_{FeatureNames[0]} .. FV_{FeatureNames[28]}
        //    Written via DiaFeatureVector.WriteTo(Span<float>)
        //
        //  Group 4 — Raw result properties for dead-feature debugging (6 cols):
        //    ObservedApexRt, LibraryRT, RtDeviationMinutes_Raw,
        //    PeakWidth_Raw, TimePointsUsed_Raw, HasPeakGroup
        //
        //  Group 5 — Phase 13 feature properties from DiaSearchResult (16 cols):
        //    MeanMassErrorPpm, MassErrorStdPpm, MaxAbsMassErrorPpm,
        //    BestFragCorrelationSum, MedianFragRefCorr, MinFragRefCorr, StdFragRefCorr,
        //    MeanSignalRatioDev, MaxSignalRatioDev, StdSignalRatioDev,
        //    SmoothedMeanFragCorr, SmoothedMinFragCorr, Log2SignalToNoise,
        //    BestFragWeightedCosine, BoundarySignalRatio, ApexToMeanRatio
        //
        //  Total: 5 + 4 + 29 + 6 + 16 = 60 columns
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