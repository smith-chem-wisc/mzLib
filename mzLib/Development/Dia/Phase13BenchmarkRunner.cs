// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
// Location: Development/Dia/Phase13BenchmarkRunner.cs

using MassSpectrometry;
using MassSpectrometry.Dia;
using Readers;
using System.Diagnostics;
using System.Globalization;

namespace Development.Dia
{
    /// <summary>
    /// Phase 13: Full 29-Feature DIA Benchmark Runner
    /// 
    /// Runs the complete DIA pipeline with Phase 13's expanded feature set (29 features):
    ///   - 3 mass accuracy features (MeanMassErrorPpm, MassErrorStdPpm, MaxAbsMassErrorPpm)
    ///   - 4 best-fragment reference curve features
    ///   - 3 signal ratio deviation features
    ///   - 2 smoothed correlation features
    ///   - 1 signal-to-noise feature
    ///   - 3 migrated features (BestFragWeightedCosine, BoundarySignalRatio, ApexToMeanRatio)
    ///   plus 13 existing features from Phase 12
    /// 
    /// Pipeline: Load → Broad Pass → RT Calibration → Calibrated Pass → (Prompts 2-4)
    /// 
    /// Runs entirely in mzLib Development project — NO MetaMorpheus dependency.
    /// </summary>
    public static class Phase13BenchmarkRunner
    {
        /// <summary>
        /// Expected feature directions for diagnostic flagging.
        /// true = higher is better for targets; false = lower is better for targets.
        /// Order matches DiaFeatureVector.FeatureNames / WriteTo().
        /// </summary>
        private static readonly bool[] HigherIsBetterForTargets = new bool[]
        {
            // [0-2]  Primary scores: ApexScore, TemporalScore, SpectralAngle
            true, true, true,
            // [3-4]  Peak correlations: PeakMeanFragCorr, PeakMinFragCorr
            true, true,
            // [5]    PeakWidth -- higher = wider elution peak = better evidence
            true,
            // [6]    CandidateCount -- lower is better (fewer ambiguous candidates)
            false,
            // [7]    LogTotalIntensity -- higher is better
            true,
            // [8]    IntensityCV -- lower is better
            false,
            // [9]    FragDetRate -- higher is better
            true,
            // [10-11] RtDeviationMinutes, RtDeviationSquared -- lower is better
            false, false,
            // [12]   TimePointsUsed -- higher is better
            true,
            // [13-15] Mass accuracy: MeanMassErrorPpm(abs), MassErrorStdPpm, MaxAbsMassErrorPpm -- lower is better
            false, false, false,
            // [16-18] BestFragCorrSum, MedianFragRefCorr, MinFragRefCorr -- higher is better
            true, true, true,
            // [19]   StdFragRefCorr -- lower is better (more consistent)
            false,
            // [20-22] MeanSigRatioDev, MaxSigRatioDev, StdSigRatioDev -- lower is better
            false, false, false,
            // [23-24] SmoothedMeanFragCorr, SmoothedMinFragCorr -- higher is better
            true, true,
            // [25]   Log2SNR -- higher is better
            true,
            // [26]   BestFragWeightedCosine -- higher is better
            true,
            // [27]   BoundarySignalRatio -- lower is better
            false,
            // [28]   ApexToMeanRatio -- higher is better
            true,
        };
        public static void RunAll(
            string rawFilePath,
            string targetMspPath,
            string decoyMspPath,
            string groundTruthTsvPath,
            string outputTsvPath)
        {
            Console.WriteLine("╔══════════════════════════════════════════════════════════════╗");
            Console.WriteLine("║    Phase 13: 29-Feature DIA Benchmark                       ║");
            Console.WriteLine("╚══════════════════════════════════════════════════════════════╝");
            Console.WriteLine();

            var totalSw = Stopwatch.StartNew();

            // ════════════════════════════════════════════════════════════
            //  Step 1: Load RT lookup
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 1: Loading RT lookup from ground truth ----------------");
            var sw = Stopwatch.StartNew();
            // ── Per-step timing capture (for Step 11 summary) ──────────
            long msRtLookup = 0, msLibLoad = 0, msRawLoad = 0, msIndexBuild = 0;
            long msBroadExtract = 0, msBroadAssemble = 0, msRtCalibration = 0;
            long msCalibratedExtract = 0, msCalibratedAssemble = 0;
            long msFeatureCompute = 0, msFdr = 0;
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
            //  Step 4: Broad pass (for RT calibration)
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 4: Broad extraction pass (RT calibration) -------------");
            var broadParams = new DiaSearchParameters
            {
                PpmTolerance = 20f,
                RtToleranceMinutes = 1.14f,   // <-- change from 5.0f
                MinFragmentsRequired = 3,
                MinScoreThreshold = 0f,
                MaxThreads = -1,
                ScoringStrategy = ScoringStrategy.TemporalCosine,
            };

            sw.Restart();
            var broadGenResult = DiaLibraryQueryGenerator.Generate(combined, index, broadParams);
            Console.WriteLine($"  Query generation: {sw.ElapsedMilliseconds}ms | {broadGenResult.Queries.Length:N0} queries " +
                $"(skipped: {broadGenResult.SkippedNoWindow} no-window, {broadGenResult.SkippedNoFragments} no-frags)");

            sw.Restart();
            ExtractionResult broadExtraction;
            using (var orchestrator = new DiaExtractionOrchestrator(index))
                broadExtraction = orchestrator.ExtractAll(broadGenResult.Queries,
                    maxDegreeOfParallelism: broadParams.EffectiveMaxThreads);
            Console.WriteLine($"  Extraction: {sw.ElapsedMilliseconds}ms | {broadExtraction.TotalDataPoints:N0} data points");

            sw.Restart();
            var broadResults = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                combined, broadGenResult, broadExtraction, broadParams, index);
            var broadAssemblyTime = sw.Elapsed;
            Console.WriteLine($"  Assembly: {broadAssemblyTime.TotalMilliseconds:F0}ms | {broadResults.Count:N0} results");

            int broadTargets = broadResults.Count(r => !r.IsDecoy);
            int broadDecoys = broadResults.Count(r => r.IsDecoy);
            Console.WriteLine($"  Targets: {broadTargets:N0} | Decoys: {broadDecoys:N0}");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 5: RT Calibration
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 5: RT Calibration --------------------------------------");

            // WORKING PATTERN from RealDataBenchmark.cs
            float anchorDpThreshold = 0.7f; // or 0.7f for cleaner anchors
            var anchors = new List<(double LibraryRt, double ObservedRt)>();

            for (int g = 0; g < broadGenResult.PrecursorGroups.Length; g++)
            {
                var group = broadGenResult.PrecursorGroups[g];
                var input = combined[group.InputIndex];

                // Skip decoys and precursors without library RT
                if (input.IsDecoy) continue;
                if (!input.RetentionTime.HasValue) continue;

                // Find the matching result to check ApexScore
                var matchingResult = broadResults.FirstOrDefault(r =>
                    r.Sequence == input.Sequence && r.ChargeState == input.ChargeState);
                if (matchingResult == null || float.IsNaN(matchingResult.ApexScore)) continue;
                if (matchingResult.ApexScore < anchorDpThreshold) continue;

                // Manually find apex RT from the raw extraction buffers
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
            Console.WriteLine($"  Anchor candidates (ApexScore >= 0.5): {anchors.Count:N0}");

            RtCalibrationModel calibration = null;
            bool useCalibrated = false;

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
                    Console.WriteLine($"  Fit time:     {fitTime.TotalMilliseconds:F1}ms");
                    Console.WriteLine($"  Model:        RT = {calibration.Slope:F4} * libRT + {calibration.Intercept:F4}");
                    Console.WriteLine($"  σ:            {calibration.SigmaMinutes:F4} min");
                    Console.WriteLine($"  R²:           {calibration.RSquared:F4}");
                    Console.WriteLine($"  Anchors used: {calibration.AnchorCount}");
                    Console.WriteLine($"  Reliable:     {(calibration.IsReliable ? "YES" : "NO")}");

                    bool usable = calibration.IsReliable || calibration.RSquared >= 0.75;

                    if (usable)
                    {
                        useCalibrated = true;
                        double halfWidth = calibration.GetMinutesWindowHalfWidth(2.0);
                        Console.WriteLine($"  Calibrated window (2σ): ±{halfWidth:F2} min " +
                            $"(vs ±{broadParams.RtToleranceMinutes:F1} min broad → " +
                            $"{broadParams.RtToleranceMinutes / halfWidth:F1}× narrower)");
                        if (!calibration.IsReliable)
                            Console.WriteLine("  NOTE: Using calibration despite !IsReliable (R² >= 0.75)");
                    }
                    else
                    {
                        Console.WriteLine("  WARNING: Calibration not reliable — falling back to broad pass results");
                    }
                }
                else
                {
                    Console.WriteLine("  WARNING: Calibration fit returned null — falling back to broad pass results");
                }
            }
            else
            {
                Console.WriteLine($"  Too few anchors ({anchors.Count}) for calibration — falling back to broad pass results");
            }
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 6: Calibrated extraction (or use broad pass)
            // ════════════════════════════════════════════════════════════
            List<DiaSearchResult> results;

            if (useCalibrated)
            {
                Console.WriteLine("--- Step 6: Calibrated extraction pass -------------------------");
                var calibratedParams = new DiaSearchParameters
                {
                    PpmTolerance = 20f,
                    RtToleranceMinutes = 5.0f,  // fallback, overridden by calibrated window
                    MinFragmentsRequired = 3,
                    ScoringStrategy = ScoringStrategy.TemporalCosine,
                    MaxThreads = -1,
                    CalibratedWindowSigmaMultiplier = 2.0,
                };

                sw.Restart();
                var calibGenResult = DiaLibraryQueryGenerator.GenerateCalibrated(
                    combined, index, calibratedParams, calibration);
                Console.WriteLine($"  Query generation: {sw.ElapsedMilliseconds}ms | {calibGenResult.Queries.Length:N0} queries");

                sw.Restart();
                ExtractionResult calibExtraction;
                using (var orchestrator = new DiaExtractionOrchestrator(index))
                    calibExtraction = orchestrator.ExtractAll(calibGenResult.Queries,
                        maxDegreeOfParallelism: calibratedParams.EffectiveMaxThreads);
                Console.WriteLine($"  Extraction: {sw.ElapsedMilliseconds}ms | {calibExtraction.TotalDataPoints:N0} data points");

                sw.Restart();
                results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    combined, calibGenResult, calibExtraction, calibratedParams, index);
                var calibAssemblyTime = sw.Elapsed;
                Console.WriteLine($"  Assembly: {calibAssemblyTime.TotalMilliseconds:F0}ms | {results.Count:N0} results");

                int calibTargets = results.Count(r => !r.IsDecoy);
                int calibDecoys = results.Count(r => r.IsDecoy);
                Console.WriteLine($"  Targets: {calibTargets:N0} | Decoys: {calibDecoys:N0}");
                Console.WriteLine($"  Δ vs broad: {results.Count - broadResults.Count:+#;-#;0} results " +
                    $"({calibExtraction.TotalDataPoints - broadExtraction.TotalDataPoints:+#;-#;0} data points)");
            }
            else
            {
                Console.WriteLine("--- Step 6: Skipped (using broad pass results) -----------------");
                results = broadResults;
            }
            Console.WriteLine();

            // ── Quick score summary ─────────────────────────────────────
            Console.WriteLine("--- Score Summary ------------------------------------------");
            var targetApex = new List<float>();
            var targetTemporal = new List<float>();
            var decoyApex = new List<float>();
            var decoyTemporal = new List<float>();

            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];
                if (r.IsDecoy)
                {
                    if (!float.IsNaN(r.ApexScore)) decoyApex.Add(r.ApexScore);
                    if (!float.IsNaN(r.TemporalScore)) decoyTemporal.Add(r.TemporalScore);
                }
                else
                {
                    if (!float.IsNaN(r.ApexScore)) targetApex.Add(r.ApexScore);
                    if (!float.IsNaN(r.TemporalScore)) targetTemporal.Add(r.TemporalScore);
                }
            }

            targetApex.Sort();
            targetTemporal.Sort();
            decoyApex.Sort();
            decoyTemporal.Sort();

            if (targetApex.Count > 0)
                Console.WriteLine($"  Target ApexScore:     median={targetApex[targetApex.Count / 2]:F4}  " +
                    $"(n={targetApex.Count:N0})");
            if (decoyApex.Count > 0)
                Console.WriteLine($"  Decoy  ApexScore:     median={decoyApex[decoyApex.Count / 2]:F4}  " +
                    $"(n={decoyApex.Count:N0})");
            if (targetTemporal.Count > 0)
                Console.WriteLine($"  Target TemporalScore: median={targetTemporal[targetTemporal.Count / 2]:F4}  " +
                    $"(n={targetTemporal.Count:N0})");
            if (decoyTemporal.Count > 0)
                Console.WriteLine($"  Decoy  TemporalScore: median={decoyTemporal[decoyTemporal.Count / 2]:F4}  " +
                    $"(n={decoyTemporal.Count:N0})");

            Console.WriteLine();
            Console.WriteLine($"  Total elapsed so far: {totalSw.Elapsed.TotalSeconds:F1}s");
            Console.WriteLine();

            // ============================================================
            //  Step 7: Compute Feature Vectors
            // ============================================================
            Console.WriteLine("--- Step 7: Computing feature vectors --------------------------");
            sw.Restart();

            if (DiaFeatureVector.ClassifierFeatureCount != 37)
            {
                Console.WriteLine($"  *** WARNING: ClassifierFeatureCount = {DiaFeatureVector.ClassifierFeatureCount}, expected 37!");
                Console.WriteLine($"  *** Action Item 5 may not be applied -- features 26-28 (migrated) may be missing.");
            }
            else
            {
                Console.WriteLine($"  ClassifierFeatureCount = {DiaFeatureVector.ClassifierFeatureCount} (confirmed: 37)");
            }

            var features = new DiaFeatureVector[results.Count];
            for (int i = 0; i < results.Count; i++)
                features[i] = DiaFeatureExtractor.ComputeFeatures(results[i], i);

            Console.WriteLine($"  Feature vectors computed: {features.Length:N0}");
            Console.WriteLine($"  Time: {sw.ElapsedMilliseconds}ms");
            Console.WriteLine();

            // ============================================================
            //  Step 8: Per-Feature Diagnostic Table
            // ============================================================
            Console.WriteLine("--- Step 8: Per-Feature Diagnostic Table -----------------------");
            sw.Restart();

            int featureCount = DiaFeatureVector.ClassifierFeatureCount;
            double[] targetSum = new double[featureCount];
            double[] decoySum = new double[featureCount];
            int[] targetCountPerFeature = new int[featureCount];
            int[] decoyCountPerFeature = new int[featureCount];
            int totalTargets = 0, totalDecoys = 0;
            Span<float> buf = stackalloc float[featureCount];

            for (int i = 0; i < results.Count; i++)
            {
                features[i].WriteTo(buf);
                bool isDecoy = results[i].IsDecoy;
                if (isDecoy) totalDecoys++; else totalTargets++;

                for (int f = 0; f < featureCount; f++)
                {
                    if (!float.IsNaN(buf[f]))
                    {
                        if (isDecoy) { decoySum[f] += buf[f]; decoyCountPerFeature[f]++; }
                        else { targetSum[f] += buf[f]; targetCountPerFeature[f]++; }
                    }
                }
            }

            Console.WriteLine($"  Total results: {results.Count:N0} (Targets: {totalTargets:N0}, Decoys: {totalDecoys:N0})");
            Console.WriteLine();

            // Table header
            Console.WriteLine("  {0,-28} | {1,11} | {2,10} | {3,10} | {4,6} | {5,6} | {6,6}",
                "Feature", "Target Mean", "Decoy Mean", "Separation", "Dir", "T NaN%", "D NaN%");
            Console.WriteLine("  " + new string('-', 96));

            // Track issues
            var allNanFeatures = new List<string>();
            var noSeparationFeatures = new List<string>();
            var reversedFeatures = new List<string>();
            var highNanFeatures = new List<string>();

            string[] featureNames = DiaFeatureVector.FeatureNames;

            for (int f = 0; f < featureCount; f++)
            {
                float targetMean = targetCountPerFeature[f] > 0
                    ? (float)(targetSum[f] / targetCountPerFeature[f]) : float.NaN;
                float decoyMean = decoyCountPerFeature[f] > 0
                    ? (float)(decoySum[f] / decoyCountPerFeature[f]) : float.NaN;

                float separation = (!float.IsNaN(targetMean) && !float.IsNaN(decoyMean))
                    ? MathF.Abs(targetMean - decoyMean) : float.NaN;

                string direction;
                if (float.IsNaN(targetMean) || float.IsNaN(decoyMean))
                    direction = "N/A";
                else if (targetMean > decoyMean)
                    direction = "T > D";
                else if (decoyMean > targetMean)
                    direction = "D > T";
                else
                    direction = "T = D";

                float targetNanPct = totalTargets > 0
                    ? 100f * (1f - (float)targetCountPerFeature[f] / totalTargets) : 100f;
                float decoyNanPct = totalDecoys > 0
                    ? 100f * (1f - (float)decoyCountPerFeature[f] / totalDecoys) : 100f;

                Console.WriteLine("  {0,-28} | {1,11} | {2,10} | {3,10} | {4,6} | {5,5:F1}% | {6,5:F1}%",
                    featureNames[f],
                    float.IsNaN(targetMean) ? "NaN" : targetMean.ToString("F4"),
                    float.IsNaN(decoyMean) ? "NaN" : decoyMean.ToString("F4"),
                    float.IsNaN(separation) ? "NaN" : separation.ToString("F4"),
                    direction,
                    targetNanPct,
                    decoyNanPct);

                string name = featureNames[f];

                if (targetCountPerFeature[f] == 0 && decoyCountPerFeature[f] == 0)
                    allNanFeatures.Add(name);

                if (!float.IsNaN(separation) && separation < 0.01f)
                    noSeparationFeatures.Add(name);

                if (!float.IsNaN(targetMean) && !float.IsNaN(decoyMean) && f < HigherIsBetterForTargets.Length)
                {
                    bool higherBetter = HigherIsBetterForTargets[f];
                    bool isReversed = higherBetter
                        ? (decoyMean > targetMean + 0.001f)
                        : (targetMean > decoyMean + 0.001f);
                    if (isReversed)
                        reversedFeatures.Add(name);
                }

                if (targetNanPct > 50f || decoyNanPct > 50f)
                    highNanFeatures.Add(name);
            }

            // Issue summary
            Console.WriteLine();
            Console.WriteLine("  -- Feature Issue Summary -----------------------------------------");

            if (allNanFeatures.Count > 0)
            {
                Console.WriteLine($"  *** ALL_NAN ({allNanFeatures.Count}) -- DATA FLOW BUG:");
                foreach (var name in allNanFeatures)
                    Console.WriteLine($"      - {name}");
            }
            else
                Console.WriteLine("  ALL_NAN: none (good -- all features have data)");

            if (noSeparationFeatures.Count > 0)
            {
                Console.WriteLine($"  *** NO_SEPARATION ({noSeparationFeatures.Count}) -- NOT DISCRIMINATIVE (|sep| < 0.01):");
                foreach (var name in noSeparationFeatures)
                    Console.WriteLine($"      - {name}");
            }
            else
                Console.WriteLine("  NO_SEPARATION: none (good -- all features show T/D separation)");

            if (reversedFeatures.Count > 0)
            {
                Console.WriteLine($"  *** REVERSED ({reversedFeatures.Count}) -- UNEXPECTED DIRECTION:");
                foreach (var name in reversedFeatures)
                    Console.WriteLine($"      - {name}");
            }
            else
                Console.WriteLine("  REVERSED: none (good -- all features in expected direction)");

            if (highNanFeatures.Count > 0)
            {
                Console.WriteLine($"  *** HIGH_NAN ({highNanFeatures.Count}) -- >50% NaN, IMPUTATION MAY HURT:");
                foreach (var name in highNanFeatures)
                    Console.WriteLine($"      - {name}");
            }
            else
                Console.WriteLine("  HIGH_NAN: none (good -- all features have >50% valid values)");

            // Category summary
            Console.WriteLine();
            Console.WriteLine("  -- Feature Category Summary --------------------------------------");
            PrintCategorySummary("Mass accuracy [13-15]", 13, 15, targetCountPerFeature, decoyCountPerFeature, targetSum, decoySum);
            PrintCategorySummary("Best-fragment [16-19]", 16, 19, targetCountPerFeature, decoyCountPerFeature, targetSum, decoySum);
            PrintCategorySummary("Signal ratio [20-22]", 20, 22, targetCountPerFeature, decoyCountPerFeature, targetSum, decoySum);
            PrintCategorySummary("Smoothed corr [23-24]", 23, 24, targetCountPerFeature, decoyCountPerFeature, targetSum, decoySum);
            PrintCategorySummary("Signal-to-noise [25]", 25, 25, targetCountPerFeature, decoyCountPerFeature, targetSum, decoySum);
            PrintCategorySummary("Migrated [26-28]", 26, 28, targetCountPerFeature, decoyCountPerFeature, targetSum, decoySum);

            Console.WriteLine();
            Console.WriteLine($"  Step 8 time: {sw.ElapsedMilliseconds}ms");
            Console.WriteLine($"  Total elapsed: {totalSw.Elapsed.TotalSeconds:F1}s");
            Console.WriteLine();

            // ═══════════════════════════════════════════════════════════════
            //  Step 9: Iterative FDR Estimation (GBT classifier)
            // ═══════════════════════════════════════════════════════════════
            sw.Restart();
            Console.WriteLine("\n═══ Step 9: FDR Estimation (GBT) ═══");

            var fdrResult = DiaFdrEngine.RunIterativeFdr(
                results, features, DiaClassifierType.LinearDiscriminant);

            sw.Stop();
            msFdr = sw.ElapsedMilliseconds;

            Console.WriteLine($"  Classifier type:   {fdrResult.ClassifierType}");
            Console.WriteLine($"  Iterations:        {fdrResult.IterationsCompleted}");
            Console.WriteLine($"  1% FDR IDs:        {fdrResult.IdentificationsAt1PctFdr:N0}");
            Console.WriteLine($"  FDR time:          {msFdr:N0} ms");

            Console.WriteLine("\n  ── Per-Iteration Diagnostics ──");
            foreach (var diag in fdrResult.Diagnostics)
            {
                DiaFdrEngine.PrintDiagnostics(diag);
            }
            // ═══════════════════════════════════════════════════════════════
            //  Step 10: FDR Summary — Q-Value Threshold Comparison
            // ═══════════════════════════════════════════════════════════════
            Console.WriteLine("\n═══ Step 10: FDR Summary ═══");

            float[] thresholds = { 0.001f, 0.005f, 0.01f, 0.05f, 0.10f };
            int[] phase12Baseline = { 5548, 9560, 11448, 19514, 27403 };

            // Count target IDs at each threshold
            int[] phase13Counts = new int[thresholds.Length];
            for (int i = 0; i < results.Count; i++)
            {
                if (results[i].IsDecoy) continue;
                var fdr = results[i].FdrInfo;
                if (fdr == null) continue;

                for (int t = 0; t < thresholds.Length; t++)
                {
                    if (fdr.QValue <= thresholds[t])
                        phase13Counts[t]++;
                }
            }

            Console.WriteLine("  Q-value threshold  |  Phase 13 IDs  |  Phase 12 IDs  |  Delta");
            Console.WriteLine("  ──────────────────────────────────────────────────────────────");
            for (int t = 0; t < thresholds.Length; t++)
            {
                int delta = phase13Counts[t] - phase12Baseline[t];
                string sign = delta >= 0 ? "+" : "";
                Console.WriteLine(
                    $"  q <= {thresholds[t]:F3}           |  {phase13Counts[t],6:N0}         |  {phase12Baseline[t],6:N0}         |  {sign}{delta:N0}");
            }

            // Target/decoy breakdown and classifier score statistics
            int totalTargets2 = 0, totalDecoys2 = 0;
            double targetScoreSum = 0, decoyScoreSum = 0;
            for (int i = 0; i < results.Count; i++)
            {
                if (results[i].IsDecoy)
                {
                    totalDecoys2++;
                    decoyScoreSum += results[i].ClassifierScore;
                }
                else
                {
                    totalTargets2++;
                    targetScoreSum += results[i].ClassifierScore;
                }
            }

            float meanTargetScore = totalTargets2 > 0 ? (float)(targetScoreSum / totalTargets2) : 0f;
            float meanDecoyScore = totalDecoys2 > 0 ? (float)(decoyScoreSum / totalDecoys2) : 0f;
            float scoreGap = meanTargetScore - meanDecoyScore;

            Console.WriteLine($"\n  Total results scored: {results.Count:N0}");
            Console.WriteLine($"  Targets: {totalTargets2:N0}   Decoys: {totalDecoys2:N0}");
            Console.WriteLine($"  Mean target classifier score: {meanTargetScore:F4}");
            Console.WriteLine($"  Mean decoy classifier score:  {meanDecoyScore:F4}");
            Console.WriteLine($"  Score gap (T - D):            {scoreGap:F4}");
            // ═══════════════════════════════════════════════════════════════
            //  Step 11: Performance Timing Summary
            // ═══════════════════════════════════════════════════════════════
            totalSw.Stop();
            Console.WriteLine("\n═══ Step 11: Timing Summary ═══");
            Console.WriteLine("  Step                      | Time");
            Console.WriteLine("  ─────────────────────────────────────");
            Console.WriteLine($"  RT lookup load            |  {msRtLookup:N0} ms");
            Console.WriteLine($"  Library load              |  {msLibLoad:N0} ms");
            Console.WriteLine($"  Raw file load + index     |  {msRawLoad:N0} ms");
            Console.WriteLine($"  Broad extraction          |  {msBroadExtract:N0} ms");
            Console.WriteLine($"  Broad assembly            |  {msBroadAssemble:N0} ms");
            Console.WriteLine($"  RT calibration fit        |  {msRtCalibration:N0} ms");
            Console.WriteLine($"  Calibrated extraction     |  {msCalibratedExtract:N0} ms");
            Console.WriteLine($"  Calibrated assembly       |  {msCalibratedAssemble:N0} ms");
            Console.WriteLine($"  Feature computation       |  {msFeatureCompute:N0} ms");
            Console.WriteLine($"  FDR estimation            |  {msFdr:N0} ms");
            Console.WriteLine($"  ─────────────────────────────────────");
            Console.WriteLine($"  TOTAL                     |  {totalSw.Elapsed.TotalSeconds:F1} s");

            // FIND (line 648):
            // Step 12 in Prompt 4

            // REPLACE WITH:
            // ═══════════════════════════════════════════════════════════════
            //  Step 12: Export Results TSV
            // ═══════════════════════════════════════════════════════════════
            if (!string.IsNullOrEmpty(outputTsvPath))
            {
                Console.WriteLine("\n═══ Step 12: TSV Export ═══");
                sw.Restart();

                ExportResultsTsv(results, features, outputTsvPath);

                sw.Stop();
                Console.WriteLine($"  Exported {results.Count:N0} rows to: {outputTsvPath}");
                Console.WriteLine($"  Export time: {sw.ElapsedMilliseconds:N0} ms");
            }
            else
            {
                Console.WriteLine("\n═══ Step 12: TSV Export (SKIPPED — no output path) ═══");
            }

            Console.WriteLine();
            Console.WriteLine("╔══════════════════════════════════════════════════════════════╗");
            Console.WriteLine("║    Phase 13 Benchmark Complete                               ║");
            Console.WriteLine("╚══════════════════════════════════════════════════════════════╝");
        }
        /// <summary>
        /// Prints a summary for a category of features (working/broken assessment).
        /// </summary>
        private static void PrintCategorySummary(
            string categoryName,
            int startIndex, int endIndex,
            int[] targetCount, int[] decoyCount,
            double[] targetSum, double[] decoySum)
        {
            int countInCategory = endIndex - startIndex + 1;
            int countWithData = 0;
            int countWithSeparation = 0;

            for (int f = startIndex; f <= endIndex; f++)
            {
                if (targetCount[f] > 0 || decoyCount[f] > 0)
                    countWithData++;

                if (targetCount[f] > 0 && decoyCount[f] > 0)
                {
                    float tMean = (float)(targetSum[f] / targetCount[f]);
                    float dMean = (float)(decoySum[f] / decoyCount[f]);
                    if (MathF.Abs(tMean - dMean) >= 0.01f)
                        countWithSeparation++;
                }
            }

            string status = countWithData == 0 ? "BROKEN (all NaN)"
                : countWithSeparation == countInCategory ? "WORKING"
                : countWithSeparation > 0 ? "PARTIAL"
                : "WEAK (no separation)";

            Console.WriteLine($"  {categoryName,-25} {countWithSeparation}/{countInCategory} discriminative -> {status}");
        }
        // ════════════════════════════════════════════════════════════════
        //  TSV Export — one row per result, all scores + Phase 13 features
        // ════════════════════════════════════════════════════════════════
        private static void ExportResultsTsv(
            List<DiaSearchResult> results,
            DiaFeatureVector[] features,
            string path)
        {
            var ic = CultureInfo.InvariantCulture;
            using var w = new StreamWriter(path);

            // ── Header ──────────────────────────────────────────────────
            // Identification
            w.Write("Sequence\tCharge\tPrecursorMz\tWindowId\tIsDecoy");
            // Scores
            w.Write("\tClassifierScore\tQValue\tPeptideQValue");
            // Primary
            w.Write("\tApexScore\tTemporalScore\tSpectralAngle\tDotProductScore");
            // Correlations
            w.Write("\tMeanFragCorr\tMinFragCorr\tPeakMeanFragCorr\tPeakMinFragCorr");
            // Fragment evidence
            w.Write("\tFragmentsDetected\tFragmentsQueried\tFragDetRate\tTimePointsUsed");
            // RT
            w.Write("\tObservedApexRt\tLibraryRT\tRtDeviationMinutes");
            // Peak group
            w.Write("\tPeakDetected\tPeakWidth\tCandidateCount");
            // Phase 13 features
            w.Write("\tMeanMassErrorPpm\tMassErrorStdPpm\tMaxAbsMassErrorPpm");
            w.Write("\tBestFragCorrelationSum\tMedianFragRefCorr\tMinFragRefCorr\tStdFragRefCorr");
            w.Write("\tMeanSignalRatioDev\tMaxSignalRatioDev\tStdSignalRatioDev");
            w.Write("\tSmoothedMeanFragCorr\tSmoothedMinFragCorr\tLog2SignalToNoise");
            w.Write("\tBestFragWeightedCosine\tBoundarySignalRatio\tApexToMeanRatio");
            // Feature vector values
            for (int j = 0; j < DiaFeatureVector.ClassifierFeatureCount; j++)
                w.Write($"\tFV_{j}");
            w.WriteLine();

            // ── Rows ────────────────────────────────────────────────────
            Span<float> buf = stackalloc float[DiaFeatureVector.ClassifierFeatureCount];

            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];
                features[i].WriteTo(buf);

                double qValue = r.FdrInfo?.QValue ?? double.NaN;
                double peptideQValue = r.FdrInfo?.PeptideQValue ?? double.NaN;

                // Identification
                w.Write(r.Sequence);
                w.Write('\t'); w.Write(r.ChargeState.ToString(ic));
                w.Write('\t'); w.Write(r.PrecursorMz.ToString("F4", ic));
                w.Write('\t'); w.Write(r.WindowId.ToString(ic));
                w.Write('\t'); w.Write(r.IsDecoy ? "True" : "False");

                // Scores
                w.Write('\t'); WriteFloat(w, r.ClassifierScore, "F6", ic);
                w.Write('\t'); WriteDouble(w, qValue, "F6", ic);
                w.Write('\t'); WriteDouble(w, peptideQValue, "F6", ic);

                // Primary
                w.Write('\t'); WriteFloat(w, r.ApexScore, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.TemporalScore, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.SpectralAngle, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.DotProductScore, "F4", ic);

                // Correlations
                w.Write('\t'); WriteFloat(w, r.MeanFragCorr, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.MinFragCorr, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.PeakMeanFragCorr, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.PeakMinFragCorr, "F4", ic);

                // Fragment evidence
                w.Write('\t'); w.Write(r.FragmentsDetected.ToString(ic));
                w.Write('\t'); w.Write(r.FragmentsQueried.ToString(ic));
                w.Write('\t'); WriteFloat(w, r.FragDetRate, "F4", ic);
                w.Write('\t'); w.Write(r.TimePointsUsed.ToString(ic));

                // RT
                w.Write('\t'); WriteFloat(w, r.ObservedApexRt, "F4", ic);
                w.Write('\t');
                if (r.LibraryRetentionTime.HasValue)
                    w.Write(r.LibraryRetentionTime.Value.ToString("F4", ic));
                else
                    w.Write("NA");
                w.Write('\t'); WriteFloat(w, r.RtDeviationMinutes, "F4", ic);

                // Peak group
                bool hasPeak = r.DetectedPeakGroup.HasValue;
                w.Write('\t'); w.Write(hasPeak ? "True" : "False");
                w.Write('\t'); WriteFloat(w, r.PeakWidth, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.CandidateCount, "F0", ic);

                // Phase 13 features
                w.Write('\t'); WriteFloat(w, r.MeanMassErrorPpm, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.MassErrorStdPpm, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.MaxAbsMassErrorPpm, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.BestFragCorrelationSum, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.MedianFragRefCorr, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.MinFragRefCorr, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.StdFragRefCorr, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.MeanSignalRatioDev, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.MaxSignalRatioDev, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.StdSignalRatioDev, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.SmoothedMeanFragCorr, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.SmoothedMinFragCorr, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.Log2SignalToNoise, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.BestFragWeightedCosine, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.BoundarySignalRatio, "F4", ic);
                w.Write('\t'); WriteFloat(w, r.ApexToMeanRatio, "F4", ic);

                // Feature vector values (FV_0 through FV_28)
                for (int j = 0; j < DiaFeatureVector.ClassifierFeatureCount; j++)
                {
                    w.Write('\t');
                    if (float.IsNaN(buf[j]))
                        w.Write("NA");
                    else
                        w.Write(buf[j].ToString("G6", ic));
                }

                w.WriteLine();
            }
        }

        private static void WriteFloat(StreamWriter w, float value, string format, CultureInfo ic)
        {
            if (float.IsNaN(value))
                w.Write("NA");
            else
                w.Write(value.ToString(format, ic));
        }

        private static void WriteDouble(StreamWriter w, double value, string format, CultureInfo ic)
        {
            if (double.IsNaN(value))
                w.Write("NA");
            else
                w.Write(value.ToString(format, ic));
        }
    }
}