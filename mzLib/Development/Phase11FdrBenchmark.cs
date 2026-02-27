// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
// Location: Development/Dia/Phase11FdrBenchmark.cs

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
    /// Phase 11: Target-Decoy FDR Estimation Benchmark
    /// 
    /// Loads target and decoy Koina .msp libraries, runs the full DIA extraction pipeline
    /// on the combined library, then runs iterative semi-supervised FDR estimation.
    /// 
    /// Reports per-iteration diagnostics and exports final results with q-values to TSV.
    /// </summary>
    public static class Phase11FdrBenchmark
    {
        public static void RunAll(
            string rawFilePath,
            string targetMspPath,
            string decoyMspPath,
            string groundTruthTsvPath,
            string outputTsvPath)
        {
            Console.WriteLine("╔══════════════════════════════════════════════════════════════╗");
            Console.WriteLine("║    Phase 11: Target-Decoy FDR Estimation Benchmark          ║");
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
            // Re-wrap with IsDecoy = true
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

            // Concatenate
            var combined = new List<LibraryPrecursorInput>(targets.Count + decoys.Count);
            combined.AddRange(targets);
            combined.AddRange(decoys);

            int hasRtTarget = 0, hasRtDecoy = 0;
            for (int i = 0; i < targets.Count; i++)
                if (targets[i].RetentionTime.HasValue) hasRtTarget++;
            for (int i = 0; i < decoys.Count; i++)
                if (decoys[i].RetentionTime.HasValue) hasRtDecoy++;

            Console.WriteLine($"  Combined: {combined.Count:N0} precursors");
            Console.WriteLine($"  Targets with RT: {hasRtTarget:N0} ({100.0 * hasRtTarget / targets.Count:F1}%)");
            Console.WriteLine($"  Decoys with RT:  {hasRtDecoy:N0} ({100.0 * hasRtDecoy / decoys.Count:F1}%)");
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
            //  Step 4: Extraction pipeline (combined target+decoy)
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
            Console.WriteLine($"  Skipped (no window): {genResult.SkippedNoWindow}");

            sw.Restart();
            using var orchestrator = new DiaExtractionOrchestrator(index);
            var extractionResult = orchestrator.ExtractAll(genResult.Queries,
                maxDegreeOfParallelism: parameters.EffectiveMaxThreads);
            Console.WriteLine($"  Extraction: {sw.ElapsedMilliseconds}ms | {extractionResult.TotalDataPoints:N0} data points");

            sw.Restart();
            var results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                combined, genResult, extractionResult, parameters);
            Console.WriteLine($"  Temporal scoring: {sw.ElapsedMilliseconds}ms | {results.Count:N0} results");

            int nTargetResults = 0, nDecoyResults = 0;
            for (int i = 0; i < results.Count; i++)
            {
                if (results[i].IsDecoy) nDecoyResults++;
                else nTargetResults++;
            }
            Console.WriteLine($"  Results: {nTargetResults:N0} targets + {nDecoyResults:N0} decoys");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 5: Compute feature vectors
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 5: Computing feature vectors --------------------------");
            sw.Restart();
            var features = new DiaFeatureVector[results.Count];
            for (int i = 0; i < results.Count; i++)
                features[i] = DiaFeatureExtractor.ComputeFeatures(results[i], i);
            Console.WriteLine($"  Feature computation: {sw.ElapsedMilliseconds}ms | {features.Length:N0} vectors");

            // Diagnostic: check RT deviation distributions for targets vs decoys
            PrintRtDeviationDiagnostics(results, features);
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 6: Iterative FDR estimation
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 6: Iterative FDR estimation ---------------------------");
            Console.WriteLine();
            sw.Restart();

            var fdrResult = DiaFdrEngine.RunIterativeFdr(
                results, features,
                maxIterations: 5,
                convergenceThreshold: 0.01f,
                idCountConvergenceThreshold: 0.01f,
                l2Lambda: 5e-3f,
                learningRate: 0.05f,
                maxEpochs: 300);

            Console.WriteLine($"  FDR estimation: {sw.ElapsedMilliseconds}ms");
            Console.WriteLine($"  Iterations: {fdrResult.IterationsCompleted}");
            Console.WriteLine($"  Final 1% FDR IDs: {fdrResult.IdentificationsAt1PctFdr:N0}");
            Console.WriteLine();

            // Per-iteration diagnostics
            Console.WriteLine("--- Per-Iteration Diagnostics ----------------------------------");
            Console.WriteLine();
            foreach (var diag in fdrResult.Diagnostics)
            {
                DiaFdrEngine.PrintDiagnostics(diag);
                Console.WriteLine();
            }

            // ════════════════════════════════════════════════════════════
            //  Step 7: Summary
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- FDR Summary ------------------------------------------------");
            Console.WriteLine();

            float[] thresholds = { 0.001f, 0.005f, 0.01f, 0.05f, 0.10f };
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

            Console.WriteLine("  Q-value threshold  |  Target IDs");
            Console.WriteLine("  ────────────────────────────────");
            for (int t = 0; t < thresholds.Length; t++)
                Console.WriteLine($"  q ≤ {thresholds[t]:F3}           |  {counts[t]:N0}");
            Console.WriteLine();

            PrintScoreDistributions(results);

            // ════════════════════════════════════════════════════════════
            //  Step 8: Export TSV
            // ════════════════════════════════════════════════════════════
            if (!string.IsNullOrEmpty(outputTsvPath))
            {
                Console.WriteLine("--- Step 8: Exporting results TSV ------------------------------");
                ExportResultsTsv(results, features, outputTsvPath);
                Console.WriteLine($"  Exported: {outputTsvPath}");
                Console.WriteLine();
            }

            totalSw.Stop();
            Console.WriteLine($"═══ Phase 11 complete. Total time: {totalSw.Elapsed.TotalSeconds:F1}s ═══");
            Console.WriteLine();
        }

        // ════════════════════════════════════════════════════════════════
        //  RT Deviation Diagnostics
        // ════════════════════════════════════════════════════════════════

        private static void PrintRtDeviationDiagnostics(
            List<DiaSearchResult> results, DiaFeatureVector[] features)
        {
            Console.WriteLine();
            Console.WriteLine("  RT deviation diagnostics (post-capping):");

            var targetRtDevs = new List<float>();
            var decoyRtDevs = new List<float>();
            int targetNoRt = 0, decoyNoRt = 0;

            for (int i = 0; i < results.Count; i++)
            {
                if (results[i].IsDecoy)
                {
                    decoyRtDevs.Add(features[i].RtDeviationMinutes);
                    if (!results[i].LibraryRetentionTime.HasValue) decoyNoRt++;
                }
                else
                {
                    targetRtDevs.Add(features[i].RtDeviationMinutes);
                    if (!results[i].LibraryRetentionTime.HasValue) targetNoRt++;
                }
            }

            targetRtDevs.Sort();
            decoyRtDevs.Sort();

            if (targetRtDevs.Count > 0)
                Console.WriteLine($"    Targets: n={targetRtDevs.Count:N0}  median={targetRtDevs[targetRtDevs.Count / 2]:F3}  " +
                    $"max={targetRtDevs[targetRtDevs.Count - 1]:F3}  noRT={targetNoRt}");
            if (decoyRtDevs.Count > 0)
                Console.WriteLine($"    Decoys:  n={decoyRtDevs.Count:N0}  median={decoyRtDevs[decoyRtDevs.Count / 2]:F3}  " +
                    $"max={decoyRtDevs[decoyRtDevs.Count - 1]:F3}  noRT={decoyNoRt}");
        }

        // ════════════════════════════════════════════════════════════════
        //  Score Distributions
        // ════════════════════════════════════════════════════════════════

        private static void PrintScoreDistributions(List<DiaSearchResult> results)
        {
            var targetScores = new List<float>();
            var decoyScores = new List<float>();

            for (int i = 0; i < results.Count; i++)
            {
                if (float.IsNaN(results[i].ClassifierScore)) continue;
                if (results[i].IsDecoy)
                    decoyScores.Add(results[i].ClassifierScore);
                else
                    targetScores.Add(results[i].ClassifierScore);
            }

            targetScores.Sort();
            decoyScores.Sort();

            Console.WriteLine("  Score distributions:");
            if (targetScores.Count > 0)
                Console.WriteLine($"    Targets (n={targetScores.Count:N0}): " +
                    $"median={targetScores[targetScores.Count / 2]:F4}  " +
                    $"Q25={targetScores[targetScores.Count / 4]:F4}  " +
                    $"Q75={targetScores[3 * targetScores.Count / 4]:F4}");
            if (decoyScores.Count > 0)
                Console.WriteLine($"    Decoys  (n={decoyScores.Count:N0}): " +
                    $"median={decoyScores[decoyScores.Count / 2]:F4}  " +
                    $"Q25={decoyScores[decoyScores.Count / 4]:F4}  " +
                    $"Q75={decoyScores[3 * decoyScores.Count / 4]:F4}");
            Console.WriteLine();
        }

        // ════════════════════════════════════════════════════════════════
        //  TSV Export
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
                w.Write('\t'); w.Write(r.ApexDotProductScore.ToString("F4"));
                w.Write('\t'); w.Write(r.TemporalCosineScore.ToString("F4"));
                w.Write('\t'); w.Write(r.SpectralAngleScore.ToString("F4"));
                w.Write('\t'); w.Write(r.MeanFragmentCorrelation.ToString("F4"));
                w.Write('\t'); w.Write(r.MinFragmentCorrelation.ToString("F4"));
                w.Write('\t'); w.Write(r.FragmentDetectionRate.ToString("F4"));
                w.Write('\t'); w.Write(r.FragmentsDetected);
                w.Write('\t'); w.Write(r.FragmentsQueried);
                w.Write('\t'); w.Write(r.TimePointsUsed);
                w.Write('\t'); w.Write(r.ObservedApexRt.ToString("F4"));
                w.Write('\t'); w.Write(r.LibraryRetentionTime.HasValue
                    ? r.LibraryRetentionTime.Value.ToString("F4") : "NA");
                w.Write('\t'); w.Write(features[i].RtDeviationMinutes.ToString("F4"));

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