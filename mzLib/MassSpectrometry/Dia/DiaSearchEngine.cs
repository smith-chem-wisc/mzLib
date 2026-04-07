// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// Location: MassSpectrometry/Dia/DiaSearchEngine.cs

using System;
using System.Collections.Generic;
using System.Diagnostics;
using MassSpectrometry.Dia.Calibration;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Production entry point for the full DIA search pipeline.
    ///
    /// Pipeline:
    ///   1. Build DiaScanIndex from raw scans
    ///   2. Bootstrap calibration → LowessRtModel (library RT → experimental RT)
    ///   3. GenerateFromLowess → per-precursor calibrated RT windows
    ///   4. ExtractAll → fragment XICs
    ///   5. AssembleResultsWithTemporalScoring → 38-feature DiaSearchResult list
    ///   6. ComputeChimericScores → feature [33]
    ///   7. ComputeFeatures → DiaFeatureVector[] (includes MS1 features [29-32])
    ///   8. RunIterativeFdr → q-values on every result
    ///
    /// DiaEngine in MetaMorpheus is a hollow shell that calls DiaSearchEngine.Run().
    /// </summary>
    public static class DiaSearchEngine
    {
        /// <summary>
        /// Result of a complete DIA search run.
        /// </summary>
        public sealed class SearchResult
        {
            /// <summary>All scored results with q-values. Targets and decoys.</summary>
            public List<DiaSearchResult> Results { get; }

            /// <summary>Feature vectors parallel to Results.</summary>
            public DiaFeatureVector[] Features { get; }

            /// <summary>The LOWESS calibration model used for this run.</summary>
            public LowessRtModel CalibrationModel { get; }

            /// <summary>Bootstrap calibration details (anchor counts, per-step results).</summary>
            public BootstrapCalibrationResult BootstrapCalibration { get; }

            /// <summary>FDR engine output (classifier, diagnostics, ID count at 1%).</summary>
            public DiaFdrEngine.FdrResult FdrResult { get; }

            /// <summary>The scan index built during this run. Caller must dispose.</summary>
            public DiaScanIndex ScanIndex { get; }

            /// <summary>Total wall-clock time for the full pipeline.</summary>
            public TimeSpan TotalTime { get; }

            internal SearchResult(
                List<DiaSearchResult> results,
                DiaFeatureVector[] features,
                LowessRtModel calibrationModel,
                BootstrapCalibrationResult bootstrapCalibration,
                DiaFdrEngine.FdrResult fdrResult,
                DiaScanIndex scanIndex,
                TimeSpan totalTime)
            {
                Results = results;
                Features = features;
                CalibrationModel = calibrationModel;
                BootstrapCalibration = bootstrapCalibration;
                FdrResult = fdrResult;
                ScanIndex = scanIndex;
                TotalTime = totalTime;
            }
        }

        /// <summary>
        /// Runs the complete DIA search pipeline.
        /// </summary>
        /// <param name="scans">
        /// All scans from the raw file — MS1 and MS2. Must not be filtered to MS2-only;
        /// MS1 scans are required for features [29-32].
        /// </param>
        /// <param name="precursors">
        /// Target + decoy library precursors. Targets must come before decoys.
        /// Each precursor must have RetentionTime set (library RT in consistent units).
        /// Precursors without RetentionTime are skipped during query generation.
        /// </param>
        /// <param name="parameters">Search parameters (ppm tolerance, threads, etc.).</param>
        /// <param name="classifierType">
        /// Classifier for iterative FDR. Default: LinearDiscriminant.
        /// LDA outperforms GBT on this feature distribution (13,629 vs 6,470 IDs at 1% FDR
        /// on HeLa benchmark). GBT score distribution collapses near zero on this data.
        /// </param>
        /// <param name="progressReporter">
        /// Optional. Receives one-line status messages. Suitable for GUI log or console.
        /// </param>
        /// <returns>
        /// A <see cref="SearchResult"/> containing scored results, feature vectors,
        /// calibration model, and FDR diagnostics.
        /// The caller is responsible for disposing <see cref="SearchResult.ScanIndex"/>.
        /// </returns>
        public static SearchResult Run(
            MsDataScan[] scans,
            IReadOnlyList<LibraryPrecursorInput> precursors,
            DiaSearchParameters parameters,
            DiaClassifierType classifierType,
            Action<string> progressReporter = null)
        {
            if (scans == null) throw new ArgumentNullException(nameof(scans));
            if (precursors == null) throw new ArgumentNullException(nameof(precursors));
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));

            void Log(string msg) => progressReporter?.Invoke(msg);

            var totalSw = Stopwatch.StartNew();

            // ── Step 1: Build scan index ─────────────────────────────────────
            Log("[DiaSearch] Building scan index...");
            var scanIndex = DiaScanIndexBuilder.Build(scans);
            Log($"[DiaSearch] Scan index: {scanIndex.Ms1ScanCount} MS1, " +
                $"{scanIndex.ScanCount} MS2, {scanIndex.WindowCount} windows, " +
                $"{scanIndex.TotalPeakCount:N0} peaks");

            if (scanIndex.ScanCount == 0)
            {
                Log("[DiaSearch] No MS2 scans found — aborting.");
                return new SearchResult(
                    new List<DiaSearchResult>(), Array.Empty<DiaFeatureVector>(),
                    null, null, default, scanIndex, totalSw.Elapsed);
            }

            // ── Step 2: Bootstrap calibration ────────────────────────────────
            Log("[DiaSearch] Running bootstrap calibration...");
            var calibSw = Stopwatch.StartNew();

            var bootstrap = DiaBootstrapCalibrator.RunBidirectional(
                precursors, scanIndex,
                ppmTolerance: parameters.PpmTolerance,
                progressReporter: progressReporter);

            calibSw.Stop();
            Log($"[DiaSearch] Calibration done in {calibSw.Elapsed.TotalSeconds:F1}s  " +
                $"anchors={bootstrap.AllAnchors.Count}  " +
                $"R²={bootstrap.LowessModel?.RSquared:F4}  " +
                $"σ={bootstrap.LowessModel?.SigmaMinutes:F4} min");

            if (bootstrap.LowessModel == null)
            {
                Log("[DiaSearch] Calibration failed (insufficient anchors) — aborting.");
                scanIndex.Dispose();
                return new SearchResult(
                    new List<DiaSearchResult>(), Array.Empty<DiaFeatureVector>(),
                    null, bootstrap, default, scanIndex, totalSw.Elapsed);
            }

            var lowess = bootstrap.LowessModel;

            // ── Step 3: Generate calibrated queries ──────────────────────────
            Log("[DiaSearch] Generating calibrated queries...");
            var genResult = DiaLibraryQueryGenerator.GenerateFromLowess(
                precursors, scanIndex, parameters.PpmTolerance, lowess);

            Log($"[DiaSearch] Queries: {genResult.Queries.Length:N0}  " +
                $"precursors: {genResult.PrecursorGroups.Length:N0}  " +
                $"skipped (no window): {genResult.SkippedNoWindow}  " +
                $"skipped (no RT/frags): {genResult.SkippedNoFragments}");

            // ── Step 4: Extract XICs ─────────────────────────────────────────
            Log("[DiaSearch] Extracting XICs...");
            var extractSw = Stopwatch.StartNew();

            ExtractionResult extractionResult;
            using (var orchestrator = new DiaExtractionOrchestrator(scanIndex))
            {
                extractionResult = orchestrator.ExtractAll(
                    genResult.Queries, parameters.EffectiveMaxThreads);
            }

            extractSw.Stop();
            Log($"[DiaSearch] Extraction done in {extractSw.Elapsed.TotalSeconds:F1}s  " +
                $"XIC points: {extractionResult.TotalDataPoints:N0}");

            // ── Step 5: Assemble results ─────────────────────────────────────
            Log("[DiaSearch] Assembling results...");
            var assembleSw = Stopwatch.StartNew();

            var results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                precursors, genResult, extractionResult, parameters, scanIndex, lowess);

            assembleSw.Stop();
            Log($"[DiaSearch] Assembly done in {assembleSw.Elapsed.TotalSeconds:F1}s  " +
                $"results: {results.Count:N0}  " +
                $"targets: {CountTargets(results)}  decoys: {CountDecoys(results)}");


            // ── Step 6: Chimeric scores ──────────────────────────────────────
            Log("[DiaSearch] Computing chimeric scores...");
            DiaLibraryQueryGenerator.ComputeChimericScores(
                precursors, results, parameters.PpmTolerance);

            // ── Step 7: Feature vectors ──────────────────────────────────────
            Log("[DiaSearch] Computing feature vectors...");
            Log($"[DiaSearch] [MS1 diag] scanIndex.Ms1ScanCount = {scanIndex.Ms1ScanCount}");

            // Build precursor lookup once — reused for both best-fragment XIC
            // extraction (MS1Ms2Correlation) and general precursor access.
            var precursorLookup = BuildPrecursorLookup(precursors);

            var features = new DiaFeatureVector[results.Count];
            int ms1DiagCount = 0; // limit trace output to first 3 non-decoy results
            for (int i = 0; i < results.Count; i++)
            {
                var result = results[i];

                // Extract best-fragment XIC for MS1/MS2 correlation [31].
                // BestFragIndex is set by DiaBestFragmentHelper during assembly.
                // If -1 (not computed), MS1Ms2Correlation stays NaN — that is fine,
                // WriteTo applies the NaN sentinel and the classifier handles it.
                float[] xicInt = Array.Empty<float>();
                float[] xicRts = Array.Empty<float>();

                var key = (result.Sequence, result.ChargeState, result.IsDecoy);
                if (result.BestFragIndex >= 0 &&
                    precursorLookup.TryGetValue(key, out int pi) &&
                    pi < precursors.Count)
                {
                    var p = precursors[pi];
                    if (result.BestFragIndex < p.FragmentCount)
                    {
                        DiaFeatureExtractor.ExtractBestFragmentXic(
                            scanIndex,
                            p.FragmentMzs[result.BestFragIndex],
                            result.WindowId,
                            result.RtWindowStart,
                            result.RtWindowEnd,
                            parameters.PpmTolerance,
                            out xicInt,
                            out xicRts);
                    }
                }

                // [MS1 diag] Before ComputeFeatures: log what is on the result and
                // what will be passed in, for the first 3 non-decoy results only.
                if (!result.IsDecoy && ms1DiagCount < 3)
                {
                    Log($"[DiaSearch] [MS1 diag] result[{i}] {result.Sequence}/{result.ChargeState}" +
                        $"  BFI={result.BestFragIndex}" +
                        $"  xicLen={xicInt.Length}" +
                        $"  PRE-ComputeFeatures: PXI={result.PrecursorXicApexIntensity:G4}" +
                        $"  IPS={result.IsotopePatternScore:G4}" +
                        $"  Corr={result.Ms1Ms2Correlation:G4}" +
                        $"  PES={result.PrecursorElutionScore:G4}");
                }

                features[i] = DiaFeatureExtractor.ComputeFeatures(
                    result, i, scanIndex,
                    bestFragXic: xicInt,
                    bestFragXicRts: xicRts);

                // [MS1 diag] After ComputeFeatures: log what ended up in the feature vector
                // and whether the result fields changed.
                if (!result.IsDecoy && ms1DiagCount < 3)
                {
                    Span<float> fvBuf = stackalloc float[DiaFeatureVector.ClassifierFeatureCount];
                    features[i].WriteTo(fvBuf);
                    Log($"[DiaSearch] [MS1 diag] result[{i}] POST-ComputeFeatures:" +
                        $"  result.PXI={result.PrecursorXicApexIntensity:G4}" +
                        $"  fv[29]={fvBuf[29]:G4}" +
                        $"  fv[30]={fvBuf[30]:G4}" +
                        $"  fv[31]={fvBuf[31]:G4}" +
                        $"  fv[32]={fvBuf[32]:G4}");
                    ms1DiagCount++;
                }
            }

            // After AssembleResultsWithTemporalScoring, before FDR
            DiaFeatureExtractor.WritePeakWidthDiagnostics(
                results,
                @"F:\DiaBenchmark\PXD005573\Runner\peakwidth_diagnostics.tsv");

            // ── Step 8: Iterative FDR ────────────────────────────────────────
            // ── Step 8: Iterative FDR ────────────────────────────────────────
            Log($"[DiaSearch] Running iterative FDR ({classifierType})...");
            var fdrSw = Stopwatch.StartNew();

            float fwhmForConsensus = lowess != null ? (float)lowess.SigmaMinutes * 0.5f : 0.069f;

            if (bootstrap.Seed != null && bootstrap.Seed.MedianPeakFwhmMinutes > 0f)
                fwhmForConsensus = bootstrap.Seed.MedianPeakFwhmMinutes;

            var fdrResult = DiaFdrEngine.RunIterativeFdr(
                results, features, classifierType,
                afterIteration1: (res, fvs) =>
                {
                    ComputeChargeStateRtConsensus(res, fvs, fwhmForConsensus);
                    Log($"[DiaSearch] ChargeStateRtConsensus injected  fwhm={fwhmForConsensus:F3} min");
                });

            fdrSw.Stop();
            Log($"[DiaSearch] FDR done in {fdrSw.Elapsed.TotalSeconds:F1}s  " +
                $"IDs at 1% FDR: {fdrResult.IdentificationsAt1PctFdr:N0}  " +
                $"iterations: {fdrResult.IterationsCompleted}");

            totalSw.Stop();
            Log($"[DiaSearch] Total time: {totalSw.Elapsed.TotalSeconds:F1}s");

            return new SearchResult(
                results, features, lowess, bootstrap, fdrResult, scanIndex, totalSw.Elapsed);
        }

        // ── Helpers ───────────────────────────────────────────────────────────

        /// <summary>
        /// Builds a lookup from (Sequence, ChargeState, IsDecoy) → index in precursors list.
        /// Used for best-fragment XIC extraction and any other per-result precursor access.
        /// </summary>
        private static Dictionary<(string, int, bool), int> BuildPrecursorLookup(
            IReadOnlyList<LibraryPrecursorInput> precursors)
        {
            var lookup = new Dictionary<(string, int, bool), int>(precursors.Count);
            for (int i = 0; i < precursors.Count; i++)
            {
                var key = (precursors[i].Sequence, precursors[i].ChargeState, precursors[i].IsDecoy);
                lookup.TryAdd(key, i);
            }
            return lookup;
        }

        private static int CountTargets(List<DiaSearchResult> results)
        {
            int n = 0;
            for (int i = 0; i < results.Count; i++)
                if (!results[i].IsDecoy) n++;
            return n;
        }

        private static int CountDecoys(List<DiaSearchResult> results)
        {
            int n = 0;
            for (int i = 0; i < results.Count; i++)
                if (results[i].IsDecoy) n++;
            return n;
        }
        /// <summary>
        /// Computes ChargeStateRtConsensus for every result.
        ///
        /// For each bare sequence, finds the result with the highest ClassifierScore
        /// across all charge states — this is the dominant charge state and its
        /// ObservedApexRt is the RT reference.
        ///
        /// Every other charge state of the same sequence gets:
        ///   ChargeStateRtConsensus = |thisApexRt - dominantApexRt| / expectedFwhmMinutes
        ///
        /// The dominant charge state itself gets 0 (it IS the reference).
        /// Singletons (only one charge state) get 0 (no information either way).
        ///
        /// Called after FDR iteration 1 has written ClassifierScore onto results,
        /// before remaining FDR iterations use the updated feature vectors.
        /// </summary>
        // Add this private static method to DiaSearchEngine:
        private static void ComputeChargeStateRtConsensus(
            List<DiaSearchResult> results,
            DiaFeatureVector[] features,
            float expectedFwhmMinutes)
        {
            // ── Pass 1: find dominant (highest ClassifierScore) apex RT per sequence ─
            var dominantApexRt = new Dictionary<(string bare, bool isDecoy), float>(results.Count / 2);
            var dominantScore = new Dictionary<(string bare, bool isDecoy), float>(results.Count / 2);

            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];
                if (float.IsNaN(r.ObservedApexRt)) continue;

                var key = (r.Sequence, r.IsDecoy);
                float score = results[i].ClassifierScore;

                if (!dominantScore.TryGetValue(key, out float best) || score > best)
                {
                    dominantScore[key] = score;
                    dominantApexRt[key] = r.ObservedApexRt;
                }
            }

            // ── Pass 2: assign deviation in units of FWHM ────────────────────────
            float invFwhm = expectedFwhmMinutes > 0f ? 1f / expectedFwhmMinutes : 1f / 0.069f;

            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];
                if (float.IsNaN(r.ObservedApexRt)) continue;

                var key = (r.Sequence, r.IsDecoy);
                if (dominantApexRt.TryGetValue(key, out float refRt))
                {
                    float dev = MathF.Abs(r.ObservedApexRt - refRt) * invFwhm;
                    r.ChargeStateRtConsensus = dev;
                    features[i].ChargeStateRtConsensus = dev;
                }
            }
        }
    }
}