// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// Location: Development/Dia/DiaSearchRunner.cs

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
    public static class DiaSearchRunner
    {
        public static void Run(
            string rawFilePath,
            string targetMspPath,
            string decoyMspPath,
            string outputTsvPath,
            string groundTruthTsvPath = null,
            DiaClassifierType classifierType = DiaClassifierType.LinearDiscriminant)
        {
            Header("DIA Search Runner — Bootstrap Calibration Pipeline");

            var totalSw = Stopwatch.StartNew();

            // ════════════════════════════════════════════════════════════
            //  Step 1: Ground truth RT lookup (optional)
            // ════════════════════════════════════════════════════════════
            Section("1. Ground truth RT lookup");
            Dictionary<string, double> rtLookup = null;
            if (!string.IsNullOrEmpty(groundTruthTsvPath) && File.Exists(groundTruthTsvPath))
            {
                rtLookup = KoinaMspParser.BuildRtLookupFromDiannTsv(groundTruthTsvPath);
                Console.WriteLine($"  Loaded {rtLookup.Count:N0} entries from {Path.GetFileName(groundTruthTsvPath)}");
            }
            else
            {
                Console.WriteLine("  No ground truth TSV — comparison skipped.");
            }

            // ════════════════════════════════════════════════════════════
            //  Step 2: Load libraries
            // ════════════════════════════════════════════════════════════
            Section("2. Loading libraries");
            if (!File.Exists(targetMspPath)) { Error($"Target MSP not found: {targetMspPath}"); return; }
            if (!File.Exists(decoyMspPath)) { Error($"Decoy MSP not found: {decoyMspPath}"); return; }

            var sw = Stopwatch.StartNew();
            var targets = KoinaMspParser.Parse(targetMspPath, rtLookup, minIntensity: 0.05f);
            Console.WriteLine($"  Targets: {targets.Count:N0} ({sw.ElapsedMilliseconds}ms)");

            sw.Restart();
            var decoysRaw = KoinaMspParser.Parse(decoyMspPath, rtLookup, minIntensity: 0.05f);
            var decoys = new List<LibraryPrecursorInput>(decoysRaw.Count);
            for (int i = 0; i < decoysRaw.Count; i++)
            {
                var d = decoysRaw[i];
                decoys.Add(new LibraryPrecursorInput(
                    d.Sequence, d.PrecursorMz, d.ChargeState,
                    d.RetentionTime, isDecoy: true,
                    d.FragmentMzs, d.FragmentIntensities, d.IrtValue));
            }
            Console.WriteLine($"  Decoys:  {decoys.Count:N0} ({sw.ElapsedMilliseconds}ms)");

            var combined = new List<LibraryPrecursorInput>(targets.Count + decoys.Count);
            combined.AddRange(targets);
            combined.AddRange(decoys);

            // Library RT distribution
            var libRts = targets.Where(p => p.RetentionTime.HasValue)
                                .Select(p => p.RetentionTime!.Value)
                                .OrderBy(x => x).ToList();
            int withRt = combined.Count(p => p.RetentionTime.HasValue);
            Console.WriteLine($"  Combined: {combined.Count:N0}  with RT: {withRt:N0} ({withRt * 100.0 / combined.Count:F1}%)");
            if (libRts.Count > 0)
            {
                Console.WriteLine($"  [DEBUG] Library RT: min={libRts.First():F3}  " +
                                  $"median={libRts[libRts.Count / 2]:F3}  " +
                                  $"max={libRts.Last():F3} min");
                Console.WriteLine($"  [DEBUG] First 5 library RTs: " +
                    string.Join(", ", libRts.Take(5).Select(r => r.ToString("F3"))));
            }

            // ════════════════════════════════════════════════════════════
            //  Step 3: Load raw file
            // ════════════════════════════════════════════════════════════
            Section("3. Loading raw file");
            Console.WriteLine($"  Path: {rawFilePath}");
            if (!File.Exists(rawFilePath)) { Error("Raw file not found."); return; }

            sw.Restart();
            MsDataFile msDataFile = Path.GetExtension(rawFilePath).ToLowerInvariant() == ".raw"
                ? (MsDataFile)new ThermoRawFileReader(rawFilePath)
                : new Mzml(rawFilePath);
            msDataFile.LoadAllStaticData();
            var scans = msDataFile.GetAllScansList().ToArray();
            sw.Stop();
            Console.WriteLine($"  {scans.Length:N0} scans in {sw.Elapsed.TotalSeconds:F2}s  " +
                              $"(MS1={scans.Count(s => s.MsnOrder == 1):N0}, " +
                              $"MS2={scans.Count(s => s.MsnOrder == 2):N0}, " +
                              $"RT {scans.First().RetentionTime:F2}–{scans.Last().RetentionTime:F2} min)");

            // ════════════════════════════════════════════════════════════
            //  Step 4: Run DiaSearchEngine
            // ════════════════════════════════════════════════════════════
            Section($"4. Running DIA search ({classifierType})");
            Console.WriteLine();

            var parameters = new DiaSearchParameters
            {
                PpmTolerance = 20f,
                MinFragmentsRequired = 3,
                MinScoreThreshold = 0f,
                MaxThreads = -1,
            };

            DiaSearchEngine.SearchResult searchResult;
            try
            {
                searchResult = DiaSearchEngine.Run(
                    scans, combined, parameters,
                    classifierType: classifierType,
                    progressReporter: msg => Console.WriteLine("  " + msg));
            }
            catch (Exception ex)
            {
                Error($"DiaSearchEngine.Run failed: {ex.Message}");
                Console.WriteLine(ex.StackTrace);
                return;
            }

            Console.WriteLine();
            var results = searchResult.Results;
            var features = searchResult.Features;

            // ════════════════════════════════════════════════════════════
            //  Step 5: Calibration diagnostics
            // ════════════════════════════════════════════════════════════
            Section("5. Calibration diagnostics");
            var bootstrap = searchResult.BootstrapCalibration;
            var lowess = searchResult.CalibrationModel;

            if (bootstrap != null)
                Console.WriteLine($"  Anchors: {bootstrap.AllAnchors.Count:N0}" +
                    $"  seed={bootstrap.Seed?.Anchors.Count ?? 0}" +
                    $"  left={bootstrap.LeftSteps.Sum(s => s.Anchors.Count)}" +
                    $"  right={bootstrap.RightSteps.Sum(s => s.Anchors.Count)}");

            if (lowess != null)
            {
                int nf = lowess.FittedLibraryRts.Length;
                Console.WriteLine($"  LOWESS: R²={lowess.RSquared:F4}  σ={lowess.SigmaMinutes:F4} min  " +
                                  $"n={lowess.AnchorCount:N0}  fitted={nf:N0}");
                if (nf > 0)
                {
                    Console.WriteLine($"  Fitted lib RT range: {lowess.FittedLibraryRts[0]:F3} – " +
                                      $"{lowess.FittedLibraryRts[nf - 1]:F3} min");
                    Console.WriteLine($"  LOWESS samples (libRT → predRT  localSig  hw):");
                    int step = Math.Max(1, nf / 10);
                    for (int i = 0; i < nf; i += step)
                    {
                        double lib = lowess.FittedLibraryRts[i];
                        double obs = lowess.ToMinutes(lib);
                        double sig = lowess.GetLocalSigma(lib);
                        double hw = Math.Min(sig, 1.5) * 2.0;
                        Console.WriteLine($"    {lib,8:F2} → {obs,7:F3}   σ={sig,6:F4}   hw={hw,6:F4}");
                    }
                }
            }
            else
            {
                Console.WriteLine("  WARN: No LOWESS model.");
            }

            // RT window diagnostics
            if (results.Count > 0)
            {
                Console.WriteLine();
                var ww = results.Select(r => r.RtWindowEnd - r.RtWindowStart).OrderBy(x => x).ToList();
                Console.WriteLine($"  [DEBUG] RT window widths — " +
                                  $"min={ww.First():F3}  p25={ww[ww.Count / 4]:F3}  " +
                                  $"median={ww[ww.Count / 2]:F3}  p75={ww[3 * ww.Count / 4]:F3}  " +
                                  $"max={ww.Last():F3} min");

                int validPeaks = results.Count(r => r.DetectedPeakGroup.HasValue && r.DetectedPeakGroup.Value.IsValid);
                Console.WriteLine($"  [DEBUG] Valid peak groups: {validPeaks:N0} / {results.Count:N0} " +
                                  $"({validPeaks * 100.0 / results.Count:F1}%)");

                var tp = results.Select(r => (float)r.TimePointsUsed).OrderBy(x => x).ToList();
                Console.WriteLine($"  [DEBUG] TimePointsUsed — " +
                                  $"min={tp.First():F0}  median={tp[tp.Count / 2]:F0}  max={tp.Last():F0}");

                Console.WriteLine($"  [DEBUG] First 5 results (libRT / window / apexRT / rtDev):");
                for (int i = 0; i < Math.Min(5, results.Count); i++)
                {
                    var r = results[i];
                    Console.WriteLine($"    [{i}] libRT={r.LibraryRetentionTime:F3}  " +
                                      $"window=[{r.RtWindowStart:F3},{r.RtWindowEnd:F3}]  " +
                                      $"apexRT={r.ObservedApexRt:F3}  " +
                                      $"rtDev={r.RtDeviationMinutes:F3}  " +
                                      $"apexScore={r.ApexScore:F4}  " +
                                      $"peakValid={r.DetectedPeakGroup?.IsValid}");
                }
            }

            // ── [DEBUG] MS1 feature pipeline trace ───────────────────────────
            // Traces every gate in the MS1 path for the first 5 non-decoy results.
            // Each line shows exactly where the path fails if MS1 features are NaN.
            Console.WriteLine();
            Console.WriteLine("  [DEBUG] MS1 feature pipeline trace (first 5 targets):");
            Console.WriteLine($"    {"Seq",-34}  {"BFI",4}  {"LookupHit",9}  {"xicLen",6}  " +
                              $"{"ms1Pts",6}  {"ApexM0",9}  {"[29]PXI",9}  {"[30]IPS",9}  " +
                              $"{"[31]Cor",9}  {"[32]PES",9}");
            Console.WriteLine("    " + new string('─', 110));

            // Rebuild precursor lookup locally for diagnostics
            var diagLookup = new Dictionary<(string, int, bool), int>(combined.Count);
            for (int i = 0; i < combined.Count; i++)
            {
                var key2 = (combined[i].Sequence, combined[i].ChargeState, combined[i].IsDecoy);
                diagLookup.TryAdd(key2, i);
            }

            int ms1DiagShown = 0;
            for (int i = 0; i < results.Count && ms1DiagShown < 5; i++)
            {
                var r = results[i];
                if (r.IsDecoy) continue;

                var key2 = (r.Sequence ?? "", r.ChargeState, r.IsDecoy);
                bool lookupHit = diagLookup.TryGetValue(key2, out int pi);

                // Reproduce what DiaSearchEngine does to get xicInt
                float[] xicInt = Array.Empty<float>();
                float[] xicRts = Array.Empty<float>();
                if (r.BestFragIndex >= 0 && lookupHit && pi < combined.Count)
                {
                    var p2 = combined[pi];
                    if (r.BestFragIndex < p2.FragmentCount)
                    {
                        DiaFeatureExtractor.ExtractBestFragmentXic(
                            searchResult.ScanIndex,
                            p2.FragmentMzs[r.BestFragIndex],
                            r.WindowId,
                            r.RtWindowStart, r.RtWindowEnd,
                            parameters.PpmTolerance,
                            out xicInt, out xicRts);
                    }
                }

                // Reproduce what Ms1XicExtractor does for isotope XIC length
                // (call ExtractIsotopeXics to see how many MS1 points are returned)
                int ms1Pts = 0;
                float apexM0 = float.NaN;
                if (searchResult.ScanIndex.Ms1ScanCount > 0)
                {
                    Ms1XicExtractor.ExtractIsotopeXics(
                        searchResult.ScanIndex,
                        (float)r.PrecursorMz,
                        Math.Max(1, r.ChargeState),
                        r.RtWindowStart, r.RtWindowEnd,
                        parameters.PpmTolerance,
                        out float[] iRts, out float[] m0Int, out _, out _);
                    ms1Pts = iRts.Length;
                    if (m0Int.Length > 0)
                    {
                        apexM0 = 0f;
                        for (int k = 0; k < m0Int.Length; k++)
                            if (m0Int[k] > apexM0) apexM0 = m0Int[k];
                    }
                }

                string pxi = float.IsNaN(r.PrecursorXicApexIntensity) ? "NaN" : r.PrecursorXicApexIntensity.ToString("F3");
                string ips = float.IsNaN(r.IsotopePatternScore) ? "NaN" : r.IsotopePatternScore.ToString("F3");
                string cor = float.IsNaN(r.Ms1Ms2Correlation) ? "NaN" : r.Ms1Ms2Correlation.ToString("F3");
                string pes = float.IsNaN(r.PrecursorElutionScore) ? "NaN" : r.PrecursorElutionScore.ToString("F3");
                string apx = float.IsNaN(apexM0) ? "NaN" : apexM0.ToString("G4");

                Console.WriteLine($"    {r.Sequence,-34}  {r.BestFragIndex,4}  {(lookupHit ? "YES" : "NO"),9}  " +
                                  $"{xicInt.Length,6}  {ms1Pts,6}  {apx,9}  {pxi,9}  {ips,9}  {cor,9}  {pes,9}");
                ms1DiagShown++;
            }

            // ── Summary counts ────────────────────────────────────────────────
            int bfiNeg = results.Count(r => !r.IsDecoy && r.BestFragIndex < 0);
            int ms1NaN = results.Count(r => !r.IsDecoy && float.IsNaN(r.PrecursorXicApexIntensity));
            int corrNaN = results.Count(r => !r.IsDecoy && float.IsNaN(r.Ms1Ms2Correlation));
            Console.WriteLine();
            Console.WriteLine($"  [DEBUG] MS1 summary (targets only, n={results.Count(r => !r.IsDecoy):N0}):");
            Console.WriteLine($"    BestFragIndex < 0:              {bfiNeg:N0}");
            Console.WriteLine($"    PrecursorXicApexIntensity=NaN:  {ms1NaN:N0}");
            Console.WriteLine($"    Ms1Ms2Correlation=NaN:          {corrNaN:N0}");
            Console.WriteLine($"    scanIndex.Ms1ScanCount:         {searchResult.ScanIndex.Ms1ScanCount:N0}");

            // ════════════════════════════════════════════════════════════
            //  Step 6: FDR results
            // ════════════════════════════════════════════════════════════
            Section("6. FDR results");
            float[] thresholds = { 0.001f, 0.005f, 0.01f, 0.05f, 0.10f };
            int[] gbtCounts = CountIdsAtThresholds(results, thresholds);

            Console.WriteLine($"  Results: {results.Count:N0}  " +
                              $"targets: {results.Count(r => !r.IsDecoy):N0}  " +
                              $"decoys: {results.Count(r => r.IsDecoy):N0}");
            Console.WriteLine();
            Console.WriteLine($"  {"Q-value",-12}  {"IDs",8}");
            Console.WriteLine("  " + new string('─', 24));
            for (int t = 0; t < thresholds.Length; t++)
                Console.WriteLine($"  q ≤ {thresholds[t]:F3}      {gbtCounts[t],8:N0}");
            Console.WriteLine();
            Console.WriteLine($"  FDR iterations: {searchResult.FdrResult.IterationsCompleted}");
            Console.WriteLine($"  IDs at 1% FDR:  {searchResult.FdrResult.IdentificationsAt1PctFdr:N0}");
            Console.WriteLine();
            Console.WriteLine("  Per-iteration diagnostics:");
            foreach (var diag in searchResult.FdrResult.Diagnostics)
                DiaFdrEngine.PrintDiagnostics(diag);

            // ════════════════════════════════════════════════════════════
            //  Step 7: Ground truth comparison
            // ════════════════════════════════════════════════════════════
            Section("7. Ground truth comparison");

            if (rtLookup != null && rtLookup.Count > 0)
            {
                int gtTotal = rtLookup.Count;
                int ourAt1Pct = gbtCounts[2];
                Console.WriteLine($"  DIA-NN (1% FDR): {gtTotal:N0}");
                Console.WriteLine($"  Ours   (1% FDR): {ourAt1Pct:N0}");
                Console.WriteLine($"  Recall:          {ourAt1Pct * 100.0 / gtTotal:F1}%");
                Console.WriteLine();

                // Build our 1% FDR set
                var ourSet = new HashSet<string>();
                for (int i = 0; i < results.Count; i++)
                {
                    var r = results[i];
                    if (r.IsDecoy || r.FdrInfo == null || r.FdrInfo.QValue > 0.01) continue;
                    ourSet.Add(r.Sequence + "/" + r.ChargeState);
                }

                int overlap = ourSet.Count(k => rtLookup.ContainsKey(k));
                int ourOnly = ourSet.Count(k => !rtLookup.ContainsKey(k));
                int gtOnly = rtLookup.Keys.Count(k => !ourSet.Contains(k));

                Console.WriteLine($"  In both:         {overlap:N0} ({overlap * 100.0 / Math.Max(ourAt1Pct, 1):F1}% of ours)");
                Console.WriteLine($"  Ours only:       {ourOnly:N0}  (false positives or novel)");
                Console.WriteLine($"  GT only (missed):{gtOnly:N0}");
                Console.WriteLine();

                // RT deviation for overlapping IDs
                var inGtDev = new List<float>();
                var notGtDev = new List<float>();
                for (int i = 0; i < results.Count; i++)
                {
                    var r = results[i];
                    if (r.IsDecoy || r.FdrInfo == null || r.FdrInfo.QValue > 0.01) continue;
                    if (float.IsNaN(r.RtDeviationMinutes)) continue;
                    string key = r.Sequence + "/" + r.ChargeState;
                    (rtLookup.ContainsKey(key) ? inGtDev : notGtDev).Add(MathF.Abs(r.RtDeviationMinutes));
                }
                if (inGtDev.Count > 0)
                {
                    inGtDev.Sort();
                    Console.WriteLine($"  |RtDev| in GT:   " +
                        $"median={inGtDev[inGtDev.Count / 2]:F3}  " +
                        $"p95={inGtDev[(int)(inGtDev.Count * 0.95)]:F3} min");
                }
                if (notGtDev.Count > 0)
                {
                    notGtDev.Sort();
                    Console.WriteLine($"  |RtDev| not GT:  " +
                        $"median={notGtDev[notGtDev.Count / 2]:F3}  " +
                        $"p95={notGtDev[(int)(notGtDev.Count * 0.95)]:F3} min");
                }
                Console.WriteLine();

                // Sample 20 missed GT IDs
                Console.WriteLine("  [DEBUG] First 20 DIA-NN IDs we missed:");
                Console.WriteLine($"    {"Key",-42}  {"GT_RT",7}  {"InOurResults",12}  {"BestScore",10}  {"BestQv",8}");
                Console.WriteLine("    " + new string('─', 85));
                int shown = 0;
                foreach (var kvp in rtLookup)
                {
                    if (ourSet.Contains(kvp.Key)) continue;
                    string inRes = "no", score = "n/a", qv = "n/a";
                    for (int i = 0; i < results.Count; i++)
                    {
                        var r = results[i];
                        if (r.IsDecoy) continue;
                        if (r.Sequence + "/" + r.ChargeState == kvp.Key)
                        {
                            inRes = "yes";
                            score = r.ClassifierScore.ToString("F4");
                            qv = r.FdrInfo?.QValue.ToString("F4") ?? "null";
                            break;
                        }
                    }
                    Console.WriteLine($"    {kvp.Key,-42}  {kvp.Value,7:F3}  {inRes,12}  {score,10}  {qv,8}");
                    if (++shown >= 20) break;
                }
                Console.WriteLine();

                // Sample 10 of our 1% IDs with RT info
                Console.WriteLine("  [DEBUG] Sample of our q≤0.01 IDs with RT info:");
                Console.WriteLine($"    {"Key",-42}  {"LibRT",7}  {"PredRT",7}  {"ApexRT",7}  {"RtDev",6}  {"InGT",5}");
                Console.WriteLine("    " + new string('─', 85));
                shown = 0;
                for (int i = 0; i < results.Count && shown < 10; i++)
                {
                    var r = results[i];
                    if (r.IsDecoy || r.FdrInfo == null || r.FdrInfo.QValue > 0.01) continue;
                    string key = r.Sequence + "/" + r.ChargeState;
                    bool inGt = rtLookup.ContainsKey(key);
                    double pred = (lowess != null && r.LibraryRetentionTime.HasValue)
                        ? lowess.ToMinutes(r.LibraryRetentionTime.Value) : double.NaN;
                    Console.WriteLine($"    {key,-42}  " +
                        $"{(r.LibraryRetentionTime.HasValue ? r.LibraryRetentionTime.Value.ToString("F3") : "null"),7}  " +
                        $"{(double.IsNaN(pred) ? "null" : pred.ToString("F3")),7}  " +
                        $"{(float.IsNaN(r.ObservedApexRt) ? "NaN" : r.ObservedApexRt.ToString("F3")),7}  " +
                        $"{(float.IsNaN(r.RtDeviationMinutes) ? "NaN" : r.RtDeviationMinutes.ToString("F3")),6}  " +
                        $"{(inGt ? "YES" : "no"),5}");
                    shown++;
                }
            }
            else
            {
                Console.WriteLine("  Skipped (no ground truth).");
                Console.WriteLine();
                Console.WriteLine("  [DEBUG] Sample q≤0.01 targets:");
                Console.WriteLine($"    {"Sequence",-42}  {"LibRT",7}  {"ApexRT",7}  {"RtDev",6}  {"Apex",6}  {"Peak",4}");
                Console.WriteLine("    " + new string('─', 82));
                int shown = 0;
                for (int i = 0; i < results.Count && shown < 15; i++)
                {
                    var r = results[i];
                    if (r.IsDecoy || r.FdrInfo == null || r.FdrInfo.QValue > 0.01) continue;
                    Console.WriteLine($"    {r.Sequence,-42}  " +
                        $"{(r.LibraryRetentionTime.HasValue ? r.LibraryRetentionTime.Value.ToString("F3") : "null"),7}  " +
                        $"{(float.IsNaN(r.ObservedApexRt) ? "NaN" : r.ObservedApexRt.ToString("F3")),7}  " +
                        $"{(float.IsNaN(r.RtDeviationMinutes) ? "NaN" : r.RtDeviationMinutes.ToString("F3")),6}  " +
                        $"{(float.IsNaN(r.ApexScore) ? "NaN" : r.ApexScore.ToString("F4")),6}  " +
                        $"{(r.DetectedPeakGroup?.IsValid.ToString() ?? "null"),4}");
                    shown++;
                }
            }

            // ════════════════════════════════════════════════════════════
            //  Step 8: Feature diagnostics
            // ════════════════════════════════════════════════════════════
            Section("8. Feature diagnostics");
            Console.WriteLine("  All targets vs all decoys:");
            PrintFeatureDiagnostics(results, features);
            Console.WriteLine();
            Console.WriteLine("  q≤0.01 targets vs all decoys:");
            PrintFeatureDiagnosticsFiltered(results, features, 0.01f);

            // ════════════════════════════════════════════════════════════
            //  Step 9: Classifier comparison
            // ════════════════════════════════════════════════════════════
            Section("9. Classifier comparison");

            Console.WriteLine("  Running LDA...");
            var ldaFdr = DiaFdrEngine.RunIterativeFdr(results, features, DiaClassifierType.LinearDiscriminant);
            int[] ldaCounts = CountIdsAtThresholds(results, thresholds);
            Console.WriteLine($"  LDA: {ldaCounts[2]:N0} IDs at 1% FDR  ({ldaFdr.IterationsCompleted} iters)");
            foreach (var diag in ldaFdr.Diagnostics) DiaFdrEngine.PrintDiagnostics(diag);

            Console.WriteLine();
            Console.WriteLine("  Running Neural Network...");
            var nnFdr = DiaFdrEngine.RunIterativeFdr(results, features, DiaClassifierType.NeuralNetwork);
            int[] nnCounts = CountIdsAtThresholds(results, thresholds);
            Console.WriteLine($"  NN:  {nnCounts[2]:N0} IDs at 1% FDR  ({nnFdr.IterationsCompleted} iters)");

            // Restore LDA scores for TSV export
            Console.WriteLine();
            Console.WriteLine("  Re-running LDA for TSV export...");
            DiaFdrEngine.RunIterativeFdr(results, features, DiaClassifierType.LinearDiscriminant);
            gbtCounts = CountIdsAtThresholds(results, thresholds);

            Console.WriteLine();
            Console.WriteLine($"  {"Classifier",-16}  {"q≤0.001",8}  {"q≤0.005",8}  {"q≤0.010",8}  {"q≤0.050",8}  {"q≤0.100",8}");
            Console.WriteLine("  " + new string('─', 68));
            PrintRow("GBT", gbtCounts);
            PrintRow("LDA", ldaCounts);
            PrintRow("NN", nnCounts);
            if (rtLookup != null)
                PrintRow("DIA-NN GT", new[] { -1, -1, rtLookup.Count, -1, -1 });

            // ════════════════════════════════════════════════════════════
            //  Step 10: Timing
            // ════════════════════════════════════════════════════════════
            totalSw.Stop();
            Section("10. Timing");
            Console.WriteLine($"  DiaSearchEngine: {searchResult.TotalTime.TotalSeconds:F1}s");
            Console.WriteLine($"  Total (incl load + 3 classifiers): {totalSw.Elapsed.TotalSeconds:F1}s");

            // ════════════════════════════════════════════════════════════
            //  Step 11: TSV export
            // ════════════════════════════════════════════════════════════
            if (!string.IsNullOrEmpty(outputTsvPath))
            {
                Section("11. TSV export");
                try
                {
                    string dir = Path.GetDirectoryName(outputTsvPath);
                    if (!string.IsNullOrEmpty(dir) && !Directory.Exists(dir))
                        Directory.CreateDirectory(dir);
                    ExportResultsTsv(results, features, outputTsvPath, "GBT");
                    Console.WriteLine($"  Written: {outputTsvPath}");
                }
                catch (Exception ex) { Error($"TSV export failed: {ex.Message}"); }
            }

            Console.WriteLine();
            Console.WriteLine(new string('═', 72));
            Console.WriteLine("  Done.");
            Console.WriteLine(new string('═', 72));
        }

        // ════════════════════════════════════════════════════════════════
        //  Feature diagnostics
        // ════════════════════════════════════════════════════════════════

        private static void PrintFeatureDiagnostics(
            List<DiaSearchResult> results, DiaFeatureVector[] features)
        {
            int n = DiaFeatureVector.ClassifierFeatureCount;
            var tSum = new double[n]; var dSum = new double[n];
            int tN = 0, dN = 0;
            Span<float> buf = stackalloc float[n];
            for (int i = 0; i < results.Count; i++)
            {
                features[i].WriteTo(buf);
                if (results[i].IsDecoy) { for (int f = 0; f < n; f++) { if (!float.IsNaN(buf[f])) dSum[f] += buf[f]; } dN++; }
                else { for (int f = 0; f < n; f++) { if (!float.IsNaN(buf[f])) tSum[f] += buf[f]; } tN++; }
            }
            Console.WriteLine($"  (n_targets={tN:N0}  n_decoys={dN:N0})");
            Console.WriteLine($"  {"Feature",-28}  {"T mean",9}  {"D mean",9}  {"Sep",7}");
            Console.WriteLine("  " + new string('─', 60));
            for (int f = 0; f < n; f++)
            {
                float tm = tN > 0 ? (float)(tSum[f] / tN) : float.NaN;
                float dm = dN > 0 ? (float)(dSum[f] / dN) : float.NaN;
                float sep = (!float.IsNaN(tm) && !float.IsNaN(dm)) ? MathF.Abs(tm - dm) : 0f;
                Console.WriteLine($"  [{f,2}] {DiaFeatureVector.FeatureNames[f],-24}  " +
                                  $"{tm,9:F4}  {dm,9:F4}  {sep,7:F4}");
            }
        }

        private static void PrintFeatureDiagnosticsFiltered(
            List<DiaSearchResult> results, DiaFeatureVector[] features, float qThresh)
        {
            int n = DiaFeatureVector.ClassifierFeatureCount;
            var tSum = new double[n]; var dSum = new double[n];
            int tN = 0, dN = 0;
            Span<float> buf = stackalloc float[n];
            for (int i = 0; i < results.Count; i++)
            {
                bool confTarget = !results[i].IsDecoy
                    && results[i].FdrInfo != null
                    && results[i].FdrInfo.QValue <= qThresh;
                if (!confTarget && !results[i].IsDecoy) continue;
                features[i].WriteTo(buf);
                if (results[i].IsDecoy) { for (int f = 0; f < n; f++) { if (!float.IsNaN(buf[f])) dSum[f] += buf[f]; } dN++; }
                else { for (int f = 0; f < n; f++) { if (!float.IsNaN(buf[f])) tSum[f] += buf[f]; } tN++; }
            }
            Console.WriteLine($"  (n_targets={tN:N0}  n_decoys={dN:N0})");
            Console.WriteLine($"  {"Feature",-28}  {"T mean",9}  {"D mean",9}  {"Sep",7}");
            Console.WriteLine("  " + new string('─', 60));
            for (int f = 0; f < n; f++)
            {
                float tm = tN > 0 ? (float)(tSum[f] / tN) : float.NaN;
                float dm = dN > 0 ? (float)(dSum[f] / dN) : float.NaN;
                float sep = (!float.IsNaN(tm) && !float.IsNaN(dm)) ? MathF.Abs(tm - dm) : 0f;
                Console.WriteLine($"  [{f,2}] {DiaFeatureVector.FeatureNames[f],-24}  " +
                                  $"{tm,9:F4}  {dm,9:F4}  {sep,7:F4}");
            }
        }

        // ════════════════════════════════════════════════════════════════
        //  TSV export
        // ════════════════════════════════════════════════════════════════

        private static void ExportResultsTsv(
            List<DiaSearchResult> results, DiaFeatureVector[] features,
            string path, string classifierType)
        {
            var inv = CultureInfo.InvariantCulture;
            using var w = new StreamWriter(path);
            w.Write("Sequence\tCharge\tPrecursorMz\tWindowId\tIsDecoy");
            w.Write("\tClassifierScore\tQValue\tClassifierType");
            for (int j = 0; j < DiaFeatureVector.ClassifierFeatureCount; j++)
                w.Write("\tFV_" + DiaFeatureVector.FeatureNames[j]);
            w.Write("\tObservedApexRt\tLibraryRT\tRtDeviationMinutes\tRtWindowStart\tRtWindowEnd");
            w.WriteLine();
            Span<float> buf = stackalloc float[DiaFeatureVector.ClassifierFeatureCount];
            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];
                double qv = r.FdrInfo?.QValue ?? 2.0;
                w.Write(r.Sequence);
                w.Write('\t'); w.Write(r.ChargeState.ToString(inv));
                w.Write('\t'); w.Write(r.PrecursorMz.ToString("F4", inv));
                w.Write('\t'); w.Write(r.WindowId.ToString(inv));
                w.Write('\t'); w.Write(r.IsDecoy ? "True" : "False");
                w.Write('\t'); w.Write(Ft(r.ClassifierScore, inv));
                w.Write('\t'); w.Write(qv.ToString("F6", inv));
                w.Write('\t'); w.Write(classifierType);
                features[i].WriteTo(buf);
                for (int j = 0; j < DiaFeatureVector.ClassifierFeatureCount; j++)
                { w.Write('\t'); w.Write(Ft(buf[j], inv)); }
                w.Write('\t'); w.Write(Ft(r.ObservedApexRt, inv));
                w.Write('\t'); w.Write(r.LibraryRetentionTime.HasValue
                    ? r.LibraryRetentionTime.Value.ToString("F4", inv) : "NA");
                w.Write('\t'); w.Write(Ft(r.RtDeviationMinutes, inv));
                w.Write('\t'); w.Write(r.RtWindowStart.ToString("F4", inv));
                w.Write('\t'); w.Write(r.RtWindowEnd.ToString("F4", inv));
                w.WriteLine();
            }
        }

        private static string Ft(float v, IFormatProvider fmt) =>
            float.IsNaN(v) ? "NA" : v.ToString("G6", fmt);

        // ════════════════════════════════════════════════════════════════
        //  Helpers
        // ════════════════════════════════════════════════════════════════

        private static void PrintRow(string name, int[] counts)
        {
            Console.Write($"  {name,-16}");
            for (int t = 0; t < counts.Length; t++)
                Console.Write(counts[t] < 0 ? $"  {"n/a",8}" : $"  {counts[t],8:N0}");
            Console.WriteLine();
        }

        private static int[] CountIdsAtThresholds(
            List<DiaSearchResult> results, float[] thresholds)
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

        private static void Header(string t)
        {
            Console.WriteLine();
            Console.WriteLine(new string('═', 72));
            Console.WriteLine("  " + t);
            Console.WriteLine(new string('═', 72));
            Console.WriteLine();
        }

        private static void Section(string t) =>
            Console.WriteLine($"\n── {t} {new string('─', Math.Max(0, 65 - t.Length))}");

        private static void Error(string m) =>
            Console.WriteLine("  [ERROR] " + m);
    }
}