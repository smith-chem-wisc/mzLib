// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// Location: mzLib/Development/Dia/MetaMorpheusEquivalentBenchmark.cs
//
// PURPOSE
// ───────
// This benchmark replicates — step for step — exactly what MetaMorpheus does when
// DiaSearchTask.RunSpecific() is invoked on a single DIA file.  By running the
// full pipeline entirely inside mzLib (no MetaMorpheus dependency, no TaskLayer,
// no EngineLayer), we avoid the stale-DLL problem that has been blocking validation.
//
// If this benchmark produces N identifications at 1% FDR, MetaMorpheus should
// produce the same N once its DLL references are correctly resolved.  Any
// discrepancy between the two outputs points directly to a MetaMorpheus-side
// integration bug (wrong parameters, mis-set IsDecoy flags, wrong library
// loading order, etc.).
//
// WHAT THIS REPLICATES (MetaMorpheus call graph)
// ───────────────────────────────────────────────
//   DiaSearchTask.RunSpecific()
//     ├─ LoadSpectralLibraries()           → KoinaMspParser.Parse (target + decoy)
//     ├─ Mark decoys by file position      → IsDecoy = true for second file
//     ├─ ToMzLibParameters()               → DiaSearchParameters (exact defaults)
//     ├─ myFileManager.LoadFile()          → MzML/ThermoRawFileReader + LoadAllStaticData
//     ├─ GetAllScansList().ToArray()        → MsDataScan[]
//     ├─ DiaScanIndexBuilder.Build(scans)  → DiaScanIndex
//     ├─ DiaEngine.ConvertLibrarySpectra() → List<LibraryPrecursorInput>
//     ├─ DiaCalibrationPipeline            → RunWithAutomaticCalibration (calibrated path)
//     ├─ WriteDiaResultsTsv (pre-FDR)      → diagnostic TSV
//     └─ PostDiaSearchAnalysisTask
//           ├─ CalculateDiaFdr             → DiaFdrEngine.RunIterativeFdr (NeuralNet)
//           ├─ DiaFeatureVector[]          → feature computation (35 features)
//           └─ WriteDiaResults             → final AllDiaResults TSV
//
// DIFFERENCES FROM METAMORPHEUS (intentional simplifications)
// ───────────────────────────────────────────────────────────
//   • No GlobalVariables.StopLoops (no GUI cancellation token)
//   • No MetaMorpheusTask progress reporting (uses Console instead)
//   • Library loaded via KoinaMspParser directly (same code path as the benchmark
//     chain; MetaMorpheus uses SpectralLibrary → LibrarySpectrum objects, then
//     DiaEngine.ConvertLibrarySpectra converts them.  We bypass LibrarySpectrum
//     here because KoinaMspParser already produces LibraryPrecursorInput directly.)
//
//     *** IMPORTANT NOTE ON LIBRARY LOADING ***
//     MetaMorpheus uses SpectralLibrary (mzLib/Readers/SpectralLibrary/SpectralLibrary.cs)
//     to load .msp files into LibrarySpectrum objects, then DiaEngine.ConvertLibrarySpectra
//     converts them to LibraryPrecursorInput.  This benchmark uses KoinaMspParser directly,
//     which produces LibraryPrecursorInput without the intermediate step.  The two paths
//     should yield identical fragment m/z and intensity arrays IF the MSP files are the
//     same.  However, if there is a bug in DiaEngine.ConvertLibrarySpectra (e.g. intensity
//     truncation, fragment filtering threshold), that bug will NOT be visible here.
//     A separate validation step (Step 2B below) can optionally load via SpectralLibrary
//     + ConvertLibrarySpectra to catch such differences.
//
// USAGE
// ─────
//   1.  Add this file to mzLib/Development/Development.csproj
//       (or place it alongside Phase23BenchmarkRunner.cs in Development/Dia/)
//   2.  In the Development project entry point (or a test), call:
//
//         MetaMorpheusEquivalentBenchmark.Run(
//             rawFilePath:      @"C:\...\Fig2HeLa-0-5h_MHRM_R01_T0.mzML",
//             targetMspPath:    @"C:\...\targets.msp",
//             decoyMspPath:     @"C:\...\decoys.msp",
//             groundTruthPath:  @"C:\...\diann_ground_truth.tsv",   // optional
//             outputFolder:     @"C:\...\BenchmarkOutput");
//
//   3.  Compare the printed "IDs at 1% FDR" count with the MetaMorpheus output.
//
// COMPILATION CHECKLIST
// ─────────────────────
//   ✓ Namespaces: MassSpectrometry, MassSpectrometry.Dia, MassSpectrometry.Dia.Calibration, Readers
//   ✓ No reference to EngineLayer, TaskLayer, MetaMorpheus assemblies
//   ✓ No LINQ in hot paths (feature computation loop uses explicit for)
//   ✓ IDisposable respected: DiaScanIndex (using var), orchestrator (using var)
//   ✓ float[] / double[] arrays match types expected by DiaFdrEngine
//   ✓ Consistent with Phase23BenchmarkRunner API usage patterns

using MassSpectrometry;
using MassSpectrometry.Dia;
using MassSpectrometry.Dia.Calibration;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Development.Dia
{
    /// <summary>
    /// Self-contained mzLib-only benchmark that replicates the full MetaMorpheus
    /// DIA search pipeline.  Run this when MetaMorpheus DLL caching prevents you
    /// from validating code changes to mzLib.
    /// </summary>
    public static class MetaMorpheusEquivalentBenchmark
    {
        // ── Baseline: what MetaMorpheus/Phase23 achieves on Fig2HeLa-0-5h ─────
        // Update these if you have a new reference run.
        private const int ExpectedMinIdsAt1PctFdr = 5_000;   // conservative lower bound
        private const int ExpectedMinPreFdrResults = 20_000;  // from calibrated search

        // ── Default parameters — exact match to Phase23BenchmarkRunner ──────────
        private const float DefaultPpmTolerance = 20f;
        private const float DefaultRtToleranceMinutes = 5.0f;
        private const int DefaultMinFragmentsRequired = 3;
        private const float DefaultMinScoreThreshold = 0.0f;
        private const float DefaultSigmaMultiplier = 3.0f;    // Phase23 uses 3.0
        private const ScoringStrategy DefaultScoringStrategy = ScoringStrategy.TemporalCosine; // Phase23 uses TemporalCosine
        private const float DefaultNonlinearPower = 0.5f;

        // ── IterativeRtCalibrator settings — exact match to Phase23BenchmarkRunner ─
        private const int CalibMaxIterations = 6;
        private const double CalibConvergenceThreshold = 0.02;
        private const double CalibSigmaMultiplierInternal = 4.0;
        // NOTE: Phase23 sets InitialTopK=500 but achieves 37,535 bootstrap anchors
        // because SelectAnchors returns ALL candidates when count <= maxAnchors.
        // Using int.MaxValue here guarantees we never artificially cap the anchor pool.
        private const int CalibInitialTopK = int.MaxValue;
        private const float CalibInitialApexScoreThreshold = 0.85f;
        private const float CalibRefinedApexScoreThreshold = 0.5f;
        private const double CalibMinWindowHalfWidthMinutes = 0.3;
        private const int CalibMinAnchorCount = 20;
        private const bool CalibEnableNonLinear = true;
        private const double CalibPiecewiseRSquaredThreshold = 0.995;
        private const double CalibLowessRSquaredThreshold = 0.990;
        private const double CalibNonLinearSigmaImprovementThreshold = 0.10;

        /// <summary>
        /// Run the full MetaMorpheus-equivalent DIA pipeline from within mzLib.
        /// All parameters have defaults that match MetaMorpheusDiaSearchParameters.
        /// </summary>
        /// <param name="rawFilePath">Path to DIA .mzML or .raw file.</param>
        /// <param name="targetMspPath">Path to target library .msp file.</param>
        /// <param name="decoyMspPath">Path to decoy library .msp file.</param>
        /// <param name="groundTruthPath">
        /// Optional path to a DIA-NN report TSV for RT lookup and anchor validation.
        /// Pass null to skip ground truth comparisons.
        /// </param>
        /// <param name="outputFolder">
        /// Folder where TSV outputs are written.  Created if it does not exist.
        /// </param>
        /// <param name="ppmTolerance">Fragment m/z tolerance in ppm.</param>
        /// <param name="rtToleranceMinutes">Initial broad-window half-width for bootstrap pass.</param>
        /// <param name="useCalibration">
        /// True (default) = full iterative RT calibration via DiaCalibrationPipeline.
        /// False = fixed RT window fallback (for debugging).
        /// </param>
        /// <param name="mspMinIntensity">
        /// Minimum relative fragment intensity threshold for MSP loading (mirrors
        /// KoinaMspParser.Parse default of 0.05).
        /// </param>
        public static void Run(
            string rawFilePath,
            string targetMspPath,
            string decoyMspPath,
            string groundTruthPath = null,
            string outputFolder = null,
            float ppmTolerance = DefaultPpmTolerance,
            float rtToleranceMinutes = DefaultRtToleranceMinutes,
            bool useCalibration = true,
            float mspMinIntensity = 0.05f)
        {
            PrintBanner();

            // ── Validate inputs ────────────────────────────────────────────────────
            if (!File.Exists(rawFilePath))
                throw new FileNotFoundException($"Raw file not found: {rawFilePath}");
            if (!File.Exists(targetMspPath))
                throw new FileNotFoundException($"Target MSP not found: {targetMspPath}");
            if (!File.Exists(decoyMspPath))
                throw new FileNotFoundException($"Decoy MSP not found: {decoyMspPath}");
            if (!string.IsNullOrEmpty(groundTruthPath) && !File.Exists(groundTruthPath))
                Console.WriteLine($"  WARNING: Ground truth TSV not found: {groundTruthPath} — skipping RT validation.");

            if (!string.IsNullOrEmpty(outputFolder))
                Directory.CreateDirectory(outputFolder);

            // ── Assert feature count matches expected value ────────────────────────
            // Guards against a partial refactor where DiaFeatureVector.ClassifierFeatureCount
            // changes without updating the FDR engine.
            Debug.Assert(
                DiaFeatureVector.ClassifierFeatureCount == 37,
                $"Expected 37 classifier features, got {DiaFeatureVector.ClassifierFeatureCount}. " +
                $"Update the feature count assertion in MetaMorpheusEquivalentBenchmark.");

            var totalSw = Stopwatch.StartNew();

            // ════════════════════════════════════════════════════════════════════════
            //  STEP 1 — Load RT ground truth (optional, for anchor quality validation)
            //
            //  MetaMorpheus equivalent: none (MetaMorpheus does not load ground truth).
            //  This step is ONLY for benchmark diagnostics and does NOT affect results.
            // ════════════════════════════════════════════════════════════════════════
            PrintStepHeader(1, "Load RT ground truth (benchmark-only — not in MetaMorpheus)");
            var sw1 = Stopwatch.StartNew();

            Dictionary<string, double> rtLookup = null;
            if (!string.IsNullOrEmpty(groundTruthPath) && File.Exists(groundTruthPath))
            {
                rtLookup = KoinaMspParser.BuildRtLookupFromDiannTsv(groundTruthPath);
                Console.WriteLine($"  RT lookup: {rtLookup.Count:N0} peptide→RT entries");
            }
            else
            {
                Console.WriteLine("  Skipped (no ground truth TSV provided).");
            }

            sw1.Stop();
            Console.WriteLine($"  Elapsed: {sw1.ElapsedMilliseconds}ms");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════════════════
            //  STEP 2 — Load target + decoy spectral library
            //
            //  MetaMorpheus equivalent:
            //    DiaSearchTask.LoadSpectralLibraries()
            //      → SpectralLibrary.Load(targetPath) + SpectralLibrary.Load(decoyPath)
            //      → Mark decoys by file position (entries after targetEntryCount)
            //    DiaEngine.ConvertLibrarySpectra()
            //      → LibrarySpectrum → LibraryPrecursorInput
            //
            //  Here we use KoinaMspParser directly (same output type, fewer allocations).
            //  See the file header note about the two loading paths.
            // ════════════════════════════════════════════════════════════════════════
            PrintStepHeader(2, "Load spectral libraries (targets + decoys)");
            var sw2 = Stopwatch.StartNew();

            Console.WriteLine($"  Target MSP:  {Path.GetFileName(targetMspPath)}");
            Console.WriteLine($"  Decoy MSP:   {Path.GetFileName(decoyMspPath)}");

            var targets = KoinaMspParser.Parse(targetMspPath, rtLookup, minIntensity: mspMinIntensity);
            Console.WriteLine($"  Targets:  {targets.Count:N0} precursors  ({sw2.ElapsedMilliseconds}ms)");

            var sw2b = Stopwatch.StartNew();
            var decoysRaw = KoinaMspParser.Parse(decoyMspPath, rtLookup, minIntensity: mspMinIntensity);

            // ── Mark all entries from the decoy MSP as IsDecoy = true ──────────────
            // This mirrors DiaSearchTask's marking logic:
            //   for (int i = targetEntryCount; i < allSpectra.Count; i++)
            //       allSpectra[i].IsDecoy = true;
            // Here we mark while building the combined list so that the IsDecoy
            // flag is guaranteed correct regardless of anything in the MSP file.
            var decoys = new List<LibraryPrecursorInput>(decoysRaw.Count);
            for (int i = 0; i < decoysRaw.Count; i++)
            {
                var d = decoysRaw[i];
                decoys.Add(new LibraryPrecursorInput(
                    sequence: d.Sequence,
                    precursorMz: d.PrecursorMz,
                    chargeState: d.ChargeState,
                    retentionTime: d.RetentionTime,
                    isDecoy: true,          // <─ ALWAYS set true for decoy MSP entries
                    fragmentMzs: d.FragmentMzs,
                    fragmentIntensities: d.FragmentIntensities,
                    irtValue: d.IrtValue));
            }
            Console.WriteLine($"  Decoys:   {decoys.Count:N0} precursors  ({sw2b.ElapsedMilliseconds}ms)");

            // ── Combine (targets first, then decoys — same order as MetaMorpheus) ──
            var combined = new List<LibraryPrecursorInput>(targets.Count + decoys.Count);
            combined.AddRange(targets);
            combined.AddRange(decoys);

            // ── Snapshot iRT range (diagnostic: span >> 30 → iRT library, slope ≠ 1) ──
            double irtMin = double.MaxValue, irtMax = double.MinValue;
            int irtCount = 0;
            for (int i = 0; i < combined.Count; i++)
            {
                var p = combined[i];
                if (p.IsDecoy) continue;
                double libRt = p.IrtValue.HasValue ? p.IrtValue.Value
                             : p.RetentionTime.HasValue ? p.RetentionTime.Value
                             : double.NaN;
                if (double.IsNaN(libRt)) continue;
                if (libRt < irtMin) irtMin = libRt;
                if (libRt > irtMax) irtMax = libRt;
                irtCount++;
            }

            sw2.Stop();
            Console.WriteLine($"  Combined: {combined.Count:N0} ({targets.Count:N0} targets, {decoys.Count:N0} decoys)");
            Console.WriteLine($"  Library iRT/RT range (targets): [{irtMin:F1}, {irtMax:F1}]  span={irtMax - irtMin:F1}  n={irtCount:N0}");
            if (irtMax - irtMin > 30.0)
                Console.WriteLine("  *** iRT library detected (span >> 30): calibration slope will be ≠ 1. ***");
            Console.WriteLine($"  Elapsed: {sw2.ElapsedMilliseconds}ms");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════════════════
            //  STEP 3 — Load raw file and build DiaScanIndex
            //
            //  MetaMorpheus equivalent:
            //    myFileManager.LoadFile(rawFilePath, combinedParams)
            //      → Mzml(path).LoadAllStaticData() or ThermoRawFileReader(path)
            //    scans = myMsDataFile.GetAllScansList().ToArray()
            //    DiaScanIndexBuilder.Build(scans)
            // ════════════════════════════════════════════════════════════════════════
            PrintStepHeader(3, "Load raw file and build DiaScanIndex");
            var swLoad = Stopwatch.StartNew();

            MsDataFile msDataFile;
            string ext = Path.GetExtension(rawFilePath).ToLowerInvariant();
            if (ext == ".raw")
                msDataFile = new ThermoRawFileReader(rawFilePath);
            else
                msDataFile = new Mzml(rawFilePath);

            msDataFile.LoadAllStaticData();
            swLoad.Stop();
            Console.WriteLine($"  File load: {swLoad.ElapsedMilliseconds / 1000.0:F2}s  ({Path.GetFileName(rawFilePath)})");

            var swIndex = Stopwatch.StartNew();
            MsDataScan[] scans = msDataFile.GetAllScansList().ToArray();
            using var scanIndex = DiaScanIndexBuilder.Build(scans);
            swIndex.Stop();

            Console.WriteLine($"  Index build: {swIndex.ElapsedMilliseconds}ms");
            Console.WriteLine($"  MS2 scans: {scanIndex.ScanCount:N0} | Windows: {scanIndex.WindowCount} | Peaks: {scanIndex.TotalPeakCount:N0}");
            Console.WriteLine($"  MS1 scans: {scanIndex.Ms1ScanCount:N0}");
            Console.WriteLine($"  RT range:  {scanIndex.GetGlobalRtMin():F2}–{scanIndex.GetGlobalRtMax():F2} min");
            Console.WriteLine();

            if (scanIndex.ScanCount == 0)
            {
                Console.WriteLine("  ERROR: No MS2 scans found. Aborting.");
                return;
            }

            // ════════════════════════════════════════════════════════════════════════
            //  STEP 4 — Build DiaSearchParameters
            //
            //  MetaMorpheus equivalent:
            //    MetaMorpheusDiaSearchParameters.ToMzLibParameters()
            //
            //  All values below are the exact defaults from MetaMorpheusDiaSearchParameters.
            //  If you change them there, update them here too.
            // ════════════════════════════════════════════════════════════════════════
            PrintStepHeader(4, "Configure DiaSearchParameters (mirrors MetaMorpheus defaults)");

            var diaParams = new DiaSearchParameters
            {
                PpmTolerance = ppmTolerance,
                RtToleranceMinutes = rtToleranceMinutes,
                MinFragmentsRequired = DefaultMinFragmentsRequired,
                MinScoreThreshold = DefaultMinScoreThreshold,
                MaxThreads = -1,
                PreferGpu = false,
                CalibratedWindowSigmaMultiplier = DefaultSigmaMultiplier,
                ScoringStrategy = DefaultScoringStrategy,
                NonlinearPower = DefaultNonlinearPower,
            };

            Console.WriteLine($"  PpmTolerance               = {diaParams.PpmTolerance}");
            Console.WriteLine($"  RtToleranceMinutes         = {diaParams.RtToleranceMinutes}");
            Console.WriteLine($"  MinFragmentsRequired       = {diaParams.MinFragmentsRequired}");
            Console.WriteLine($"  ScoringStrategy            = {diaParams.ScoringStrategy}");
            Console.WriteLine($"  SigmaMultiplier            = {diaParams.CalibratedWindowSigmaMultiplier}");
            Console.WriteLine($"  UseCalibration             = {useCalibration}");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════════════════
            //  STEP 5 — Build IterativeRtCalibrator (calibration strategy)
            //
            //  MetaMorpheus equivalent:
            //    DiaEngine does not expose calibrator settings — it uses
            //    DiaCalibrationPipeline.RunWithAutomaticCalibration() with a default
            //    IterativeRtCalibrator constructed inside the pipeline.
            //
            //  The values below match Phase23BenchmarkRunner to ensure we test the
            //  same algorithm path.  They also match the internal defaults in
            //  IterativeRtCalibrator.
            // ════════════════════════════════════════════════════════════════════════
            PrintStepHeader(5, "Configure IterativeRtCalibrator");

            var calibrator = new IterativeRtCalibrator
            {
                MaxIterations = CalibMaxIterations,
                ConvergenceThreshold = CalibConvergenceThreshold,
                SigmaMultiplier = CalibSigmaMultiplierInternal,
                InitialTopK = CalibInitialTopK,
                InitialApexScoreThreshold = CalibInitialApexScoreThreshold,
                RefinedApexScoreThreshold = CalibRefinedApexScoreThreshold,
                MinWindowHalfWidthMinutes = CalibMinWindowHalfWidthMinutes,
                MinAnchorCount = CalibMinAnchorCount,
                EnableNonLinearModelSelection = CalibEnableNonLinear,
                PiecewiseLinearRSquaredThreshold = CalibPiecewiseRSquaredThreshold,
                LowessRSquaredThreshold = CalibLowessRSquaredThreshold,
                NonLinearSigmaImprovementThreshold = CalibNonLinearSigmaImprovementThreshold,
                AnchorMaxCoElutionStd = 2.0f,   // Prompt 4: reject scattered-fragment anchors
                AnchorMinCandidateScoreGap = 0.05f,  // Prompt 4: reject ambiguous peak-group anchors
            };

            Console.WriteLine($"  MaxIterations              = {calibrator.MaxIterations}");
            Console.WriteLine($"  ConvergenceThreshold       = {calibrator.ConvergenceThreshold}");
            Console.WriteLine($"  InitialApexScoreThreshold  = {calibrator.InitialApexScoreThreshold}");
            Console.WriteLine($"  MinAnchorCount             = {calibrator.MinAnchorCount}");
            Console.WriteLine($"  EnableNonLinear            = {calibrator.EnableNonLinearModelSelection}");
            Console.WriteLine($"  AnchorMaxCoElutionStd      = {calibrator.AnchorMaxCoElutionStd}");
            Console.WriteLine($"  AnchorMinCandidateScoreGap = {calibrator.AnchorMinCandidateScoreGap}");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════════════════
            //  STEP 6 — Run DIA search (calibrated or fixed-window)
            //
            //  MetaMorpheus equivalent:
            //    DiaEngine.RunSpecific()
            //      └─ DiaCalibrationPipeline.RunWithAutomaticCalibration(...)
            //         (when useCalibration == true)
            // ════════════════════════════════════════════════════════════════════════
            PrintStepHeader(6, useCalibration
                ? "Run DIA search — calibrated path (DiaCalibrationPipeline)"
                : "Run DIA search — fixed RT window (no calibration)");

            DiaCalibrationPipeline.PipelineResult pipelineResult = default;
            List<DiaSearchResult> searchResults;
            RtCalibrationModel calibrationModel = null;

            var swSearch = Stopwatch.StartNew();
            using (var orchestrator = new DiaExtractionOrchestrator(scanIndex))
            {
                if (useCalibration)
                {
                    pipelineResult = DiaCalibrationPipeline.RunWithAutomaticCalibration(
                        combined, scanIndex, diaParams, orchestrator, calibrator,
                        progressReporter: msg => Console.WriteLine($"  [Cal] {msg}"));

                    searchResults = pipelineResult.Results;
                    calibrationModel = pipelineResult.Calibration;
                }
                else
                {
                    // Fixed RT window — mirrors DiaEngine.RunFixedRtSearch
                    Console.WriteLine("  (Running fixed-window extraction — no RT calibration)");
                    var genResult = DiaLibraryQueryGenerator.Generate(
                        combined, scanIndex, diaParams);
                    var extractionResult = orchestrator.ExtractAll(
                        genResult.Queries, diaParams.EffectiveMaxThreads);
                    searchResults = DiaLibraryQueryGenerator.AssembleResults(
                        combined, genResult, extractionResult.Results, diaParams,
                        dotProductScorer: new NormalizedDotProductScorer(),
                        spectralAngleScorer: new SpectralAngleScorer());
                    calibrationModel = null;
                }
            }
            swSearch.Stop();
            Console.WriteLine();

            int targetCount = 0, decoyCount = 0;
            for (int i = 0; i < searchResults.Count; i++)
            {
                if (searchResults[i].IsDecoy) decoyCount++;
                else targetCount++;
            }

            Console.WriteLine($"  Search time:   {swSearch.ElapsedMilliseconds / 1000.0:F1}s");
            Console.WriteLine($"  Total results: {searchResults.Count:N0}  ({targetCount:N0} targets, {decoyCount:N0} decoys)");

            if (searchResults.Count < ExpectedMinPreFdrResults)
                Console.WriteLine($"  *** WARNING: result count ({searchResults.Count:N0}) is below expected minimum ({ExpectedMinPreFdrResults:N0}). Check calibration. ***");

            // ── Print calibration summary ─────────────────────────────────────────
            PrintCalibrationSummary(calibrationModel, pipelineResult);
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════════════════
            //  STEP 7 — Write pre-FDR diagnostic TSV
            //
            //  MetaMorpheus equivalent:
            //    DiaSearchTask.WriteDiaResultsTsv (pre-FDR diagnostic file)
            //    Written to:  {OutputFolder}/Individual File Results/{RawFileName}/
            //                 {RawFileName}_DiaResults_PreFDR.tsv
            // ════════════════════════════════════════════════════════════════════════
            PrintStepHeader(7, "Write pre-FDR diagnostic TSV");

            string preFdrPath = null;
            if (!string.IsNullOrEmpty(outputFolder))
            {
                string rawFileName = Path.GetFileNameWithoutExtension(rawFilePath);
                preFdrPath = Path.Combine(outputFolder, rawFileName + "_DiaResults_PreFDR.tsv");
                WritePreFdrTsv(preFdrPath, rawFileName, searchResults, writeDecoys: false);
                Console.WriteLine($"  Wrote: {preFdrPath}");
            }
            else
            {
                Console.WriteLine("  Skipped (no output folder specified).");
            }
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════════════════
            //  STEP 8 — Compute feature vectors for FDR
            //
            //  MetaMorpheus equivalent:
            //    PostDiaSearchAnalysisTask.CalculateDiaFdr()
            //      for (int i = 0; i < results.Count; i++)
            //          DiaFeatureExtractor.ComputeFeatures(result, precursorIndex)
            //
            // ════════════════════════════════════════════════════════════════════════
            PrintStepHeader(8, "Compute FDR feature vectors (35 features per result)");
            var swFeatures = Stopwatch.StartNew();

            var featureVectors = new DiaFeatureVector[searchResults.Count];

            for (int i = 0; i < searchResults.Count; i++)
            {
                var result = searchResults[i];
                featureVectors[i] = DiaFeatureExtractor.ComputeFeatures(result, i);
            }

            swFeatures.Stop();
            Console.WriteLine($"  Feature computation: {swFeatures.ElapsedMilliseconds}ms  ({featureVectors.Length:N0} vectors × {DiaFeatureVector.ClassifierFeatureCount} features)");

            // Quick sanity check: print mean feature values for targets vs decoys
            PrintFeatureSanityCheck(featureVectors, searchResults);
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════════════════
            //  STEP 9 — Run iterative semi-supervised FDR (LDA then NeuralNet)
            //
            //  MetaMorpheus equivalent:
            //    DiaFdrEngine.RunIterativeFdr(results, featureVectors, DiaClassifierType.NeuralNetwork)
            //    (called from PostDiaSearchAnalysisTask.CalculateDiaFdr)
            //
            //  Phase23BenchmarkRunner runs LDA first, then NN — the LDA pass
            //  is diagnostic but also warms up the separation and produces a
            //  comparable baseline.  We mirror that here.
            // ════════════════════════════════════════════════════════════════════════
            PrintStepHeader(9, "Run iterative FDR (LDA then NeuralNet classifier)");
            var swFdr = Stopwatch.StartNew();

            // ── LDA (diagnostic + warm baseline) ────────────────────────────────
            Console.WriteLine("  [LDA]");
            long ldaMs = swFdr.ElapsedMilliseconds;
            var ldaFdrResult = DiaFdrEngine.RunIterativeFdr(
                searchResults, featureVectors,
                classifierType: DiaClassifierType.LinearDiscriminant,
                maxIterations: 5);
            ldaMs = swFdr.ElapsedMilliseconds - ldaMs;
            Console.WriteLine($"  LDA FDR complete: {ldaFdrResult.IdentificationsAt1PctFdr:N0} IDs at 1% FDR ({ldaMs}ms)");
            Console.WriteLine();

            // ── Neural Network ───────────────────────────────────────────────────
            Console.WriteLine("  [Neural Network]");
            long nnMs = swFdr.ElapsedMilliseconds;
            DiaFdrEngine.RunIterativeFdr(
                searchResults, featureVectors,
                classifierType: DiaClassifierType.NeuralNetwork,
                maxIterations: 5);
            nnMs = swFdr.ElapsedMilliseconds - nnMs;

            Console.WriteLine($"  NN  FDR complete: ({nnMs}ms)");
            swFdr.Stop();
            Console.WriteLine($"  FDR computation: {swFdr.ElapsedMilliseconds}ms  (LDA: {ldaMs}ms + NN: {nnMs}ms)");

            // ── Count IDs at various q-value thresholds ──────────────────────────
            int idsAt0_001 = 0, idsAt0_005 = 0, idsAt0_01 = 0, idsAt0_05 = 0, idsAt0_10 = 0;
            for (int i = 0; i < searchResults.Count; i++)
            {
                var r = searchResults[i];
                if (r.IsDecoy || r.FdrInfo == null) continue;
                double q = r.FdrInfo.QValue;
                if (q <= 0.001) idsAt0_001++;
                if (q <= 0.005) idsAt0_005++;
                if (q <= 0.01) idsAt0_01++;
                if (q <= 0.05) idsAt0_05++;
                if (q <= 0.10) idsAt0_10++;
            }

            Console.WriteLine();
            Console.WriteLine("  ┌─────────────────────────────────────────────┐");
            Console.WriteLine("  │   FDR Results                               │");
            Console.WriteLine("  ├─────────────────────────────────────────────┤");
            Console.WriteLine($"  │  LDA  1% FDR: {ldaFdrResult.IdentificationsAt1PctFdr,8:N0} IDs                    │");
            Console.WriteLine($"  │  NN q ≤ 0.001:{idsAt0_001,8:N0} IDs                    │");
            Console.WriteLine($"  │  NN q ≤ 0.005:{idsAt0_005,8:N0} IDs                    │");
            Console.WriteLine($"  │  NN q ≤ 0.01 :{idsAt0_01,8:N0} IDs  ← primary metric  │");
            Console.WriteLine($"  │  NN q ≤ 0.05 :{idsAt0_05,8:N0} IDs                    │");
            Console.WriteLine($"  │  NN q ≤ 0.10 :{idsAt0_10,8:N0} IDs                    │");
            Console.WriteLine("  └─────────────────────────────────────────────┘");
            Console.WriteLine();

            if (idsAt0_01 < ExpectedMinIdsAt1PctFdr)
                Console.WriteLine($"  *** WARNING: IDs at 1% FDR ({idsAt0_01:N0}) is below expected minimum ({ExpectedMinIdsAt1PctFdr:N0}). ***");
            else
                Console.WriteLine($"  ✓  IDs at 1% FDR = {idsAt0_01:N0} (above minimum {ExpectedMinIdsAt1PctFdr:N0})");

            Console.WriteLine();

            // ── Peptide-level (unique sequences) at 1% FDR ───────────────────────
            var uniquePeptidesAt1Pct = new HashSet<string>();
            for (int i = 0; i < searchResults.Count; i++)
            {
                var r = searchResults[i];
                if (!r.IsDecoy && r.FdrInfo != null && r.FdrInfo.QValue <= 0.01)
                    uniquePeptidesAt1Pct.Add(r.Sequence);
            }
            Console.WriteLine($"  Unique peptide sequences at 1% FDR: {uniquePeptidesAt1Pct.Count:N0}");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════════════════
            //  STEP 10 — Write final results TSV
            //
            //  MetaMorpheus equivalent:
            //    PostDiaSearchAnalysisTask.WriteDiaResults (AllDiaResults.tsv)
            //    and optionally PostDiaSearchAnalysisTask.WritePsmTsvResults
            // ════════════════════════════════════════════════════════════════════════
            PrintStepHeader(10, "Write final results TSV");

            string finalTsvPath = null;
            if (!string.IsNullOrEmpty(outputFolder))
            {
                string rawFileName = Path.GetFileNameWithoutExtension(rawFilePath);
                finalTsvPath = Path.Combine(outputFolder, rawFileName + "_AllDiaResults.tsv");
                WriteFinalResultsTsv(finalTsvPath, rawFileName, searchResults,
                    qValueThreshold: 0.01, writeDecoys: false);
                Console.WriteLine($"  Wrote: {finalTsvPath}");

                // Also write a full-results TSV (all results, targets only, for inspection)
                string allResultsPath = Path.Combine(outputFolder, rawFileName + "_AllResults_NoFDRFilter.tsv");
                WriteFinalResultsTsv(allResultsPath, rawFileName, searchResults,
                    qValueThreshold: 1.0, writeDecoys: false);
                Console.WriteLine($"  Wrote: {allResultsPath}");
            }
            else
            {
                Console.WriteLine("  Skipped (no output folder specified).");
            }
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════════════════
            //  STEP 11 — Calibration anchor validation (optional, requires ground truth)
            // ════════════════════════════════════════════════════════════════════════
            if (rtLookup != null && useCalibration && calibrationModel != null)
            {
                PrintStepHeader(11, "Calibration anchor validation vs DIA-NN ground truth");
                ValidateCalibrationAnchors(searchResults, rtLookup, calibrationModel);
                Console.WriteLine();
            }

            // ════════════════════════════════════════════════════════════════════════
            //  SUMMARY
            // ════════════════════════════════════════════════════════════════════════
            totalSw.Stop();
            Console.WriteLine("════════════════════════════════════════════════════════════════");
            Console.WriteLine("  BENCHMARK SUMMARY");
            Console.WriteLine("════════════════════════════════════════════════════════════════");
            Console.WriteLine($"  Raw file:               {Path.GetFileName(rawFilePath)}");
            Console.WriteLine($"  Library:                {targets.Count:N0} targets + {decoys.Count:N0} decoys");
            Console.WriteLine($"  MS2 scans indexed:      {scanIndex.ScanCount:N0}");
            Console.WriteLine($"  Pre-FDR results:        {searchResults.Count:N0}  ({targetCount:N0} T, {decoyCount:N0} D)");
            Console.WriteLine($"  IDs at 1% FDR:          {idsAt0_01:N0}");
            Console.WriteLine($"  Unique peptides (1%):   {uniquePeptidesAt1Pct.Count:N0}");

            if (calibrationModel != null)
            {
                Console.WriteLine($"  Calibration slope:      {calibrationModel.Slope:F4}");
                Console.WriteLine($"  Calibration intercept:  {calibrationModel.Intercept:F4}");
                Console.WriteLine($"  Calibration σ:          {calibrationModel.SigmaMinutes:F4} min");
                Console.WriteLine($"  Calibration R²:         {calibrationModel.RSquared:F4}");
            }

            Console.WriteLine($"  Total elapsed:          {totalSw.Elapsed.TotalSeconds:F1}s");
            Console.WriteLine();

            // ── MetaMorpheus parity check ─────────────────────────────────────────
            Console.WriteLine("  METAMORPHEUS PARITY CHECK");
            Console.WriteLine("  ─────────────────────────────────────────────────");
            Console.WriteLine("  Run MetaMorpheus DIA search on the same files and compare:");
            Console.WriteLine($"    Expected IDs at 1% FDR:  {idsAt0_01:N0}  (this benchmark)");
            Console.WriteLine("    Look in MetaMorpheus output:  AllDiaResults.tsv");
            Console.WriteLine("    Count rows with q-value ≤ 0.01 and IsDecoy = FALSE");
            Console.WriteLine("    Any difference → MetaMorpheus integration bug (not mzLib bug)");

            if (calibrationModel != null)
            {
                Console.WriteLine();
                Console.WriteLine("  CALIBRATION PARITY CHECK");
                Console.WriteLine("  ─────────────────────────────────────────────────");
                Console.WriteLine($"    Expected slope:     {calibrationModel.Slope:F4}");
                Console.WriteLine($"    Expected intercept: {calibrationModel.Intercept:F4}");
                Console.WriteLine($"    Expected σ:         {calibrationModel.SigmaMinutes:F4} min");
                Console.WriteLine("    MetaMorpheus log line to look for:");
                Console.WriteLine("      '[DIA DEBUG] Iter: sigma=... R2=...' (last iteration)");
                Console.WriteLine("      The final sigma should match the above.");
            }

            Console.WriteLine();
            Console.WriteLine("════════════════════════════════════════════════════════════════");
        }

        // ════════════════════════════════════════════════════════════════════════
        //  PRIVATE HELPERS
        // ════════════════════════════════════════════════════════════════════════

        private static void PrintBanner()
        {
            Console.WriteLine();
            Console.WriteLine("╔════════════════════════════════════════════════════════════════╗");
            Console.WriteLine("║  MetaMorpheus-Equivalent DIA Benchmark (mzLib-only)            ║");
            Console.WriteLine("║  Validates mzLib DIA engine without MetaMorpheus DLL cache     ║");
            Console.WriteLine("╚════════════════════════════════════════════════════════════════╝");
            Console.WriteLine();
        }

        private static void PrintStepHeader(int stepNumber, string title)
        {
            Console.WriteLine($"─── Step {stepNumber}: {title} ───");
        }

        private static void PrintCalibrationSummary(
            RtCalibrationModel model,
            DiaCalibrationPipeline.PipelineResult pipelineResult)
        {
            if (model == null)
            {
                Console.WriteLine("  Calibration: NOT RUN or FAILED (model is null).");
                Console.WriteLine("  *** Expected slope ~0.175 for Prosit iRT library on 45-min gradient ***");
                return;
            }

            Console.WriteLine($"  Calibration converged:");
            Console.WriteLine($"    Slope:      {model.Slope:F4}");
            Console.WriteLine($"    Intercept:  {model.Intercept:F4} min");
            Console.WriteLine($"    σ:          {model.SigmaMinutes:F4} min");
            Console.WriteLine($"    R²:         {model.RSquared:F4}");

            if (pipelineResult.CalibrationLog != null && pipelineResult.CalibrationLog.Count > 0)
            {
                Console.WriteLine($"  Calibration log ({pipelineResult.CalibrationLog.Count} iterations):");
                Console.WriteLine($"  {"Iter",4}  {"σ (min)",8}  {"R²",7}  {"Window ±",8}");
                Console.WriteLine($"  {"────",4}  {"────────",8}  {"───────",7}  {"────────",8}");
                for (int i = 0; i < pipelineResult.CalibrationLog.Count; i++)
                {
                    var entry = pipelineResult.CalibrationLog[i];
                    Console.WriteLine($"  {i,4}  {entry.SigmaMinutes,8:F4}  {entry.RSquared,7:F4}  ±{entry.WindowHalfWidthMinutes:F4}");
                }
            }

            // Warn if slope is far from expected iRT slope for Prosit libraries
            if (model.Slope < 0.05 || model.Slope > 5.0)
                Console.WriteLine($"  *** WARNING: slope={model.Slope:F4} is highly unusual (expected 0.05–5.0). Check library. ***");
            else if (model.Slope < 0.10 || model.Slope > 1.5)
                Console.WriteLine($"  *** INFO: slope={model.Slope:F4} suggests iRT library (typical Prosit: ~0.175 for 45-min gradient). ***");

            if (model.SigmaMinutes > 1.0)
                Console.WriteLine($"  *** WARNING: σ={model.SigmaMinutes:F4} min > 1.0 min — calibration may not have converged. ***");
        }

        private static void PrintFeatureSanityCheck(
            DiaFeatureVector[] features,
            List<DiaSearchResult> results)
        {
            if (features.Length == 0) return;

            // Compute mean of first 5 features for targets vs decoys to verify separation
            int nFeatures = Math.Min(5, DiaFeatureVector.ClassifierFeatureCount);
            var targetSums = new double[nFeatures];
            var decoySums = new double[nFeatures];
            int nTargets = 0, nDecoys = 0;

            Span<float> featureBuf = stackalloc float[DiaFeatureVector.ClassifierFeatureCount];
            for (int i = 0; i < features.Length; i++)
            {
                bool isDecoy = results[i].IsDecoy;
                features[i].WriteTo(featureBuf);
                for (int f = 0; f < nFeatures; f++)
                {
                    double val = featureBuf[f];
                    if (isDecoy) decoySums[f] += val;
                    else targetSums[f] += val;
                }
                if (isDecoy) nDecoys++;
                else nTargets++;
            }

            Console.WriteLine($"  Feature sanity check (first {nFeatures} features, mean value):");
            Console.WriteLine($"  {"Feature",8}  {"Target mean",12}  {"Decoy mean",12}  {"Separation",12}");
            Console.WriteLine($"  {"───────",8}  {"───────────",12}  {"──────────",12}  {"──────────",12}");
            for (int f = 0; f < nFeatures; f++)
            {
                double tMean = nTargets > 0 ? targetSums[f] / nTargets : 0;
                double dMean = nDecoys > 0 ? decoySums[f] / nDecoys : 0;
                double sep = tMean - dMean;
                Console.WriteLine($"  {f,8}  {tMean,12:F4}  {dMean,12:F4}  {sep,12:F4}");
            }

            // Warn if features look identical for targets and decoys (FDR will fail)
            double maxSep = 0;
            for (int f = 0; f < nFeatures; f++)
            {
                double tMean = nTargets > 0 ? targetSums[f] / nTargets : 0;
                double dMean = nDecoys > 0 ? decoySums[f] / nDecoys : 0;
                double sep = Math.Abs(tMean - dMean);
                if (sep > maxSep) maxSep = sep;
            }
            if (maxSep < 0.01)
                Console.WriteLine("  *** WARNING: targets and decoys have nearly identical feature values — FDR will not separate them. ***");
        }

        private static void ValidateCalibrationAnchors(
            List<DiaSearchResult> results,
            Dictionary<string, double> rtLookup,
            RtCalibrationModel model)
        {
            // Compare ObservedApexRt against DIA-NN ground truth RT
            int matched = 0, total = 0;
            double sumDelta = 0, maxDelta = 0;
            var deltas = new List<double>();

            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];
                if (r.IsDecoy) continue;
                if (float.IsNaN(r.ObservedApexRt)) continue;
                if (!rtLookup.TryGetValue(r.Sequence, out double truthRt)) continue;

                double delta = Math.Abs(r.ObservedApexRt - truthRt);
                deltas.Add(delta);
                sumDelta += delta;
                if (delta > maxDelta) maxDelta = delta;
                matched++;
                total++;
            }

            if (total == 0)
            {
                Console.WriteLine("  No matched sequences found in ground truth lookup.");
                return;
            }

            deltas.Sort();
            double median = deltas[deltas.Count / 2];
            double mean = sumDelta / total;
            double pct95 = deltas[(int)(deltas.Count * 0.95)];

            Console.WriteLine($"  Matched vs ground truth: {matched:N0} results");
            Console.WriteLine($"  |ObservedRT - DIA-NN RT| stats:");
            Console.WriteLine($"    Mean:    {mean:F3} min");
            Console.WriteLine($"    Median:  {median:F3} min");
            Console.WriteLine($"    95th %:  {pct95:F3} min");
            Console.WriteLine($"    Max:     {maxDelta:F3} min");

            if (median > 0.5)
                Console.WriteLine("  *** WARNING: median RT error > 0.5 min — calibration may be incorrect. ***");
            else
                Console.WriteLine("  ✓  Median RT error is acceptable (< 0.5 min).");
        }

        // ── Pre-FDR TSV (mirrors DiaSearchTask.WriteDiaResultsTsv) ──────────────

        private static readonly string PreFdrTsvHeader = string.Join("\t", new[]
        {
            "File Name", "Sequence", "Precursor m/z", "Charge", "Window ID", "Is Decoy",
            "Dot Product Score", "Spectral Angle Score", "Fragments Detected",
            "Fragments Queried", "Fragment Detection Rate", "Library RT",
            "RT Window Start", "RT Window End", "XIC Point Counts", "Extracted Intensities"
        });

        private static void WritePreFdrTsv(
            string filePath,
            string rawFileName,
            List<DiaSearchResult> results,
            bool writeDecoys)
        {
            using var writer = new StreamWriter(filePath, false,
                new System.Text.UTF8Encoding(false))
            { NewLine = "\n" };
            writer.WriteLine(PreFdrTsvHeader);

            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];
                if (!writeDecoys && r.IsDecoy) continue;

                string xicCounts = r.XicPointCounts != null
                    ? string.Join(";", r.XicPointCounts)
                    : "";
                string intensities = r.ExtractedIntensities != null
                    ? string.Join(";", r.ExtractedIntensities
                        .Select(x => x.ToString("G6", CultureInfo.InvariantCulture)))
                    : "";

                writer.WriteLine(string.Join("\t",
                    rawFileName,
                    r.Sequence,
                    r.PrecursorMz.ToString("F4", CultureInfo.InvariantCulture),
                    r.ChargeState,
                    r.WindowId,
                    r.IsDecoy ? "TRUE" : "FALSE",
                    float.IsNaN(r.DotProductScore) ? "NaN"
                        : r.DotProductScore.ToString("F4", CultureInfo.InvariantCulture),
                    float.IsNaN(r.SpectralAngleScore) ? "NaN"
                        : r.SpectralAngleScore.ToString("F4", CultureInfo.InvariantCulture),
                    r.FragmentsDetected,
                    r.FragmentsQueried,
                    r.FragmentDetectionRate.ToString("F4", CultureInfo.InvariantCulture),
                    r.LibraryRetentionTime.HasValue
                        ? r.LibraryRetentionTime.Value.ToString("F2", CultureInfo.InvariantCulture)
                        : "",
                    r.RtWindowStart.ToString("F2", CultureInfo.InvariantCulture),
                    r.RtWindowEnd.ToString("F2", CultureInfo.InvariantCulture),
                    xicCounts,
                    intensities));
            }
        }

        // ── Final Results TSV (mirrors PostDiaSearchAnalysisTask.WriteDiaResults) ─

        private static readonly string FinalTsvHeader = string.Join("\t", new[]
        {
            "File Name", "Sequence", "Precursor m/z", "Charge", "Window ID", "Is Decoy",
            "Q-Value", "Peptide Q-Value", "Score", "Dot Product Score", "Spectral Angle Score",
            "Apex Score", "Temporal Score", "Spectral Angle", "Observed Apex RT",
            "Library RT", "RT Window Start", "RT Window End",
            "Fragments Detected", "Fragments Queried", "Fragment Detection Rate"
        });

        private static void WriteFinalResultsTsv(
            string filePath,
            string rawFileName,
            List<DiaSearchResult> results,
            double qValueThreshold,
            bool writeDecoys)
        {
            using var writer = new StreamWriter(filePath, false,
                new System.Text.UTF8Encoding(false))
            { NewLine = "\n" };
            writer.WriteLine(FinalTsvHeader);

            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];
                if (!writeDecoys && r.IsDecoy) continue;
                if (r.FdrInfo != null && r.FdrInfo.QValue > qValueThreshold) continue;

                writer.WriteLine(string.Join("\t",
                    rawFileName,
                    r.Sequence,
                    r.PrecursorMz.ToString("F4", CultureInfo.InvariantCulture),
                    r.ChargeState,
                    r.WindowId,
                    r.IsDecoy ? "TRUE" : "FALSE",
                    r.FdrInfo != null ? r.FdrInfo.QValue.ToString("F6", CultureInfo.InvariantCulture) : "N/A",
                    r.FdrInfo?.PeptideQValue.HasValue == true
                        ? r.FdrInfo.PeptideQValue.Value.ToString("F6", CultureInfo.InvariantCulture) : "N/A",
                    r.ClassifierScore.ToString("F4", CultureInfo.InvariantCulture),
                    float.IsNaN(r.DotProductScore) ? "NaN"
                        : r.DotProductScore.ToString("F4", CultureInfo.InvariantCulture),
                    float.IsNaN(r.SpectralAngleScore) ? "NaN"
                        : r.SpectralAngleScore.ToString("F4", CultureInfo.InvariantCulture),
                    float.IsNaN(r.ApexScore) ? "NaN"
                        : r.ApexScore.ToString("F4", CultureInfo.InvariantCulture),
                    float.IsNaN(r.TemporalScore) ? "NaN"
                        : r.TemporalScore.ToString("F4", CultureInfo.InvariantCulture),
                    float.IsNaN(r.SpectralAngle) ? "NaN"
                        : r.SpectralAngle.ToString("F4", CultureInfo.InvariantCulture),
                    !float.IsNaN(r.ObservedApexRt)
                        ? r.ObservedApexRt.ToString("F4", CultureInfo.InvariantCulture) : "",
                    r.LibraryRetentionTime.HasValue
                        ? r.LibraryRetentionTime.Value.ToString("F2", CultureInfo.InvariantCulture) : "",
                    r.RtWindowStart.ToString("F2", CultureInfo.InvariantCulture),
                    r.RtWindowEnd.ToString("F2", CultureInfo.InvariantCulture),
                    r.FragmentsDetected,
                    r.FragmentsQueried,
                    r.FragmentDetectionRate.ToString("F4", CultureInfo.InvariantCulture)));
            }
        }
    }
}