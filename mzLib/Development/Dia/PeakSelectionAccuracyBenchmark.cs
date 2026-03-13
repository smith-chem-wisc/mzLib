// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
// Location: Development/Dia/PeakSelectionAccuracyBenchmark.cs

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
    /// Prompt 1 of 5 — Peak Selection Accuracy Benchmark
    ///
    /// Diagnostic tool that quantifies exactly how bad (or good) peak group
    /// selection is as a function of extraction window width.
    ///
    /// This benchmark is NOT part of the production pipeline. It is a
    /// before/after measurement harness that should be run:
    ///   - Before Prompt 2: to establish the current-detector baseline
    ///   - After Prompt 2: to measure improvement from the new detector
    ///
    /// Key question answered: At each extraction window width, what fraction
    /// of precursors have their peak group apex correctly identified?
    ///
    /// Correct  = |observedApexRt - trueRt| &lt; 0.5 min
    /// Marginal = 0.5 – 1.5 min
    /// Wrong    = > 1.5 min, or NaN
    /// </summary>
    public static class PeakSelectionAccuracyBenchmark
    {
        // ─── Classification thresholds (minutes) ────────────────────────
        private const float CorrectThreshold = 0.5f;
        private const float MarginalThreshold = 1.5f;

        // ─── Window half-widths to sweep ────────────────────────────────
        private static readonly float[] HalfWidths = { 1.0f, 2.0f, 3.0f, 5.0f };

        // ─── Benchmark precursor count targets ──────────────────────────
        private const int DecilesCount = 10;
        private const int MaxPerDecile = 20;
        private const int TotalCap = 200; // soft cap

        // ─── Classification codes ────────────────────────────────────────
        private const string Correct = "Correct";
        private const string Marginal = "Marginal";
        private const string Wrong = "Wrong";

        // ────────────────────────────────────────────────────────────────
        //  Entry point
        // ────────────────────────────────────────────────────────────────

        /// <summary>
        /// Run the peak selection accuracy benchmark across four RT window widths.
        /// </summary>
        /// <param name="rawFilePath">Path to .raw or .mzml file</param>
        /// <param name="targetMspPath">Path to target MSP library (Koina/NIST format)</param>
        /// <param name="groundTruthTsvPath">DIA-NN TSV with true per-peptide RTs</param>
        /// <param name="outputTsvPath">Where to write per-precursor detail TSV</param>
        public static void Run(
            string rawFilePath,
            string targetMspPath,
            string groundTruthTsvPath,
            string outputTsvPath)
        {
            Console.WriteLine("╔══════════════════════════════════════════════════════════════════╗");
            Console.WriteLine("║        Peak Selection Accuracy Benchmark  (Prompt 1 of 5)        ║");
            Console.WriteLine("╚══════════════════════════════════════════════════════════════════╝");
            Console.WriteLine();

            var totalSw = Stopwatch.StartNew();

            // ════════════════════════════════════════════════════════════
            //  Step 1 — Load ground truth RT lookup
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 1: Loading RT lookup from ground truth ----------------");
            var sw = Stopwatch.StartNew();

            if (string.IsNullOrEmpty(groundTruthTsvPath) || !File.Exists(groundTruthTsvPath))
                throw new FileNotFoundException($"Ground truth TSV not found: {groundTruthTsvPath}");

            Dictionary<string, double> trueRtBySequence =
                KoinaMspParser.BuildRtLookupFromDiannTsv(groundTruthTsvPath);

            Console.WriteLine($"  RT lookup entries: {trueRtBySequence.Count:N0}  ({sw.ElapsedMilliseconds} ms)");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 2 — Load target library (no decoys needed here)
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 2: Loading target library -----------------------------");
            sw.Restart();

            var allTargets = KoinaMspParser.Parse(targetMspPath, trueRtBySequence, minIntensity: 0.05f);
            Console.WriteLine($"  Targets loaded: {allTargets.Count:N0}  ({sw.ElapsedMilliseconds} ms)");

            // Diagnose key-format mismatch between library and ground truth
            int matched = 0;
            for (int i = 0; i < allTargets.Count; i++)
                if (trueRtBySequence.ContainsKey(allTargets[i].Sequence + "/" + allTargets[i].ChargeState)) matched++;

            double matchFrac = allTargets.Count > 0 ? (double)matched / allTargets.Count : 0.0;
            Console.WriteLine($"  Ground truth coverage: {matched:N0}/{allTargets.Count:N0} = {matchFrac:P1}");

            if (matchFrac < 0.5)
            {
                Console.WriteLine();
                Console.WriteLine("  *** WARNING: fewer than 50% of library sequences matched ground truth keys.");
                Console.WriteLine("  *** This likely indicates a sequence format mismatch.");
                Console.WriteLine("  *** Showing 3 example keys from each side for diagnosis:");
                int shown = 0;
                foreach (var t in allTargets)
                {
                    if (shown >= 3) break;
                    Console.WriteLine($"      Library  key: '{t.Sequence}'");
                    shown++;
                }
                shown = 0;
                foreach (var kvp in trueRtBySequence)
                {
                    if (shown >= 3) break;
                    Console.WriteLine($"      GT TSV   key: '{kvp.Key}'");
                    shown++;
                }
                Console.WriteLine();
            }
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 3 — Select up to 200 benchmark precursors
            //           Stratified across 10 iRT deciles, ≤20 per decile
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 3: Selecting benchmark precursors (stratified by iRT) -");

            // Filter to those that have a ground truth RT entry
            var withGt = new List<LibraryPrecursorInput>(allTargets.Count);
            for (int i = 0; i < allTargets.Count; i++)
                if (trueRtBySequence.ContainsKey(allTargets[i].Sequence + "/" + allTargets[i].ChargeState))
                    withGt.Add(allTargets[i]);

            Console.WriteLine($"  Precursors with ground truth RT: {withGt.Count:N0}");

            var benchmarkPrecursors = StratifyByIrtDecile(withGt, DecilesCount, MaxPerDecile, TotalCap);

            // Determine how many iRT deciles are represented
            int decilesRepresented = CountDecilesRepresented(withGt, benchmarkPrecursors, DecilesCount);

            Console.WriteLine($"  Selected {benchmarkPrecursors.Count} benchmark precursors across {decilesRepresented} iRT deciles");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 4 — Load raw file and build index (once, reused across windows)
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 4: Loading raw file and building scan index -----------");
            sw.Restart();

            MsDataFile msDataFile;
            string ext = Path.GetExtension(rawFilePath).ToLowerInvariant();
            if (ext == ".raw")
                msDataFile = new ThermoRawFileReader(rawFilePath);
            else
                msDataFile = new Mzml(rawFilePath);
            msDataFile.LoadAllStaticData();
            long loadMs = sw.ElapsedMilliseconds;

            sw.Restart();
            var scans = msDataFile.GetAllScansList().ToArray();
            using var index = DiaScanIndexBuilder.Build(scans);
            long indexMs = sw.ElapsedMilliseconds;

            Console.WriteLine($"  File load: {loadMs / 1000.0:F1}s | Index build: {indexMs}ms");
            Console.WriteLine($"  Scans: {index.ScanCount:N0} | Windows: {index.WindowCount} | Peaks: {index.TotalPeakCount:N0}");
            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 5 — Run extraction at 4 window widths and collect records
            // ════════════════════════════════════════════════════════════
            Console.WriteLine("--- Step 5: Running extraction across 4 window widths ----------");

            // All per-precursor records across all windows, for TSV export
            var allRecords = new List<PeakAccuracyRecord>(benchmarkPrecursors.Count * HalfWidths.Length);

            // Per-window summary stats (indexed by HalfWidths position)
            var summaries = new WindowSummary[HalfWidths.Length];

            for (int wi = 0; wi < HalfWidths.Length; wi++)
            {
                float halfWidth = HalfWidths[wi];
                Console.WriteLine($"  Window ±{halfWidth:F1} min ...");
                sw.Restart();

                var parameters = new DiaSearchParameters
                {
                    PpmTolerance = 20f,
                    RtToleranceMinutes = halfWidth,
                    MinFragmentsRequired = 3,
                    MinScoreThreshold = 0f,
                    MaxThreads = -1,
                    ScoringStrategy = ScoringStrategy.TemporalCosine,
                };

                // Generate queries for benchmark precursors only (not full library)
                var genResult = DiaLibraryQueryGenerator.Generate(benchmarkPrecursors, index, parameters);

                using var orchestrator = new DiaExtractionOrchestrator(index);
                var extractionResult = orchestrator.ExtractAll(
                    genResult.Queries,
                    maxDegreeOfParallelism: parameters.EffectiveMaxThreads);

                var results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                    benchmarkPrecursors, genResult, extractionResult, parameters, index);

                long elapsedMs = sw.ElapsedMilliseconds;

                // Build a quick lookup: sequence → result (take first match if dupes)
                var resultBySeq = new Dictionary<string, DiaSearchResult>(results.Count, StringComparer.Ordinal);
                for (int r = 0; r < results.Count; r++)
                {
                    string key = ResultKey(results[r]);
                    if (!resultBySeq.ContainsKey(key))
                        resultBySeq[key] = results[r];
                }

                // Collect records for every benchmark precursor
                int correct = 0, marginal = 0, wrong = 0;
                var absErrors = new List<float>(benchmarkPrecursors.Count);

                for (int p = 0; p < benchmarkPrecursors.Count; p++)
                {
                    var prec = benchmarkPrecursors[p];
                    string precKey = PrecursorKey(prec);

                    double trueRt = trueRtBySequence[prec.Sequence + "/" + prec.ChargeState];
                    double iRt = prec.IrtValue ?? prec.RetentionTime ?? double.NaN;

                    float observedApexRt = float.NaN;
                    float apexCosine = float.NaN;
                    int candidateCount = 0;
                    string classification;

                    if (resultBySeq.TryGetValue(precKey, out var res))
                    {
                        observedApexRt = res.ObservedApexRt;
                        apexCosine = res.ApexScore;
                        candidateCount = (int)Math.Round(res.CandidateCount);
                    }

                    float absError;
                    if (float.IsNaN(observedApexRt))
                    {
                        absError = 99.0f;
                        classification = Wrong;
                        wrong++;
                    }
                    else
                    {
                        absError = MathF.Abs(observedApexRt - (float)trueRt);
                        if (absError < CorrectThreshold)
                        {
                            classification = Correct;
                            correct++;
                            absErrors.Add(absError);
                        }
                        else if (absError < MarginalThreshold)
                        {
                            classification = Marginal;
                            marginal++;
                            absErrors.Add(absError);
                        }
                        else
                        {
                            classification = Wrong;
                            wrong++;
                            absErrors.Add(absError);
                        }
                    }

                    allRecords.Add(new PeakAccuracyRecord
                    {
                        Sequence = prec.Sequence,
                        ChargeState = prec.ChargeState,
                        IRt = (float)iRt,
                        TrueRt = (float)trueRt,
                        WindowHalf = halfWidth,
                        SelectedApexRt = observedApexRt,
                        AbsError = absError,
                        Classification = classification,
                        ApexCosine = apexCosine,
                        CandidateCount = candidateCount,
                    });
                }

                float medianError = ComputeMedian(absErrors);
                int total = benchmarkPrecursors.Count;

                summaries[wi] = new WindowSummary
                {
                    HalfWidth = halfWidth,
                    Total = total,
                    Correct = correct,
                    Marginal = marginal,
                    Wrong = wrong,
                    MedianError = medianError,
                };

                Console.WriteLine($"    Correct={correct}/{total}  Marginal={marginal}/{total}  " +
                                  $"Wrong={wrong}/{total}  Median err={medianError:F3} min  ({elapsedMs} ms)");
            }

            Console.WriteLine();

            // ════════════════════════════════════════════════════════════
            //  Step 6 — Print accuracy table
            // ════════════════════════════════════════════════════════════
            PrintAccuracyTable(summaries);

            // ════════════════════════════════════════════════════════════
            //  Step 7 — Print 10 worst misses at ±5 min
            // ════════════════════════════════════════════════════════════
            PrintWorstMisses(allRecords, worstWindowHalf: 5.0f, topN: 10);

            // ════════════════════════════════════════════════════════════
            //  Step 8 — Write detail TSV
            // ════════════════════════════════════════════════════════════
            WriteDetailTsv(allRecords, outputTsvPath);
            Console.WriteLine($"Detail TSV written: {outputTsvPath}");
            Console.WriteLine();

            Console.WriteLine($"Total benchmark time: {totalSw.Elapsed.TotalSeconds:F1}s");
        }

        // ────────────────────────────────────────────────────────────────
        //  Stratified selection by iRT decile
        // ────────────────────────────────────────────────────────────────

        private static List<LibraryPrecursorInput> StratifyByIrtDecile(
            List<LibraryPrecursorInput> candidates,
            int deciles,
            int maxPerDecile,
            int totalCap)
        {
            if (candidates.Count == 0)
                return new List<LibraryPrecursorInput>();

            // Sort by iRT (fall back to RetentionTime if iRT unavailable)
            var sorted = new List<LibraryPrecursorInput>(candidates);
            sorted.Sort((a, b) =>
            {
                double aKey = a.IrtValue ?? a.RetentionTime ?? 0.0;
                double bKey = b.IrtValue ?? b.RetentionTime ?? 0.0;
                return aKey.CompareTo(bKey);
            });

            int n = sorted.Count;
            var selected = new List<LibraryPrecursorInput>(Math.Min(totalCap, deciles * maxPerDecile));

            for (int d = 0; d < deciles; d++)
            {
                int lo = (int)((long)d * n / deciles);
                int hi = (int)((long)(d + 1) * n / deciles);
                int take = Math.Min(maxPerDecile, hi - lo);

                // Take evenly spaced within the decile slice to maximize coverage
                for (int k = 0; k < take; k++)
                {
                    int idx = lo + (int)((long)k * (hi - lo) / take);
                    selected.Add(sorted[idx]);
                    if (selected.Count >= totalCap) break;
                }
                if (selected.Count >= totalCap) break;
            }

            return selected;
        }

        private static int CountDecilesRepresented(
            List<LibraryPrecursorInput> withGt,
            List<LibraryPrecursorInput> selected,
            int deciles)
        {
            if (withGt.Count == 0 || selected.Count == 0) return 0;

            // Build a set of selected sequences for fast lookup
            var selSet = new HashSet<string>(selected.Count, StringComparer.Ordinal);
            foreach (var s in selected) selSet.Add(s.Sequence);

            // Sort withGt by iRT to assign decile ids
            var sorted = new List<LibraryPrecursorInput>(withGt);
            sorted.Sort((a, b) =>
            {
                double aKey = a.IrtValue ?? a.RetentionTime ?? 0.0;
                double bKey = b.IrtValue ?? b.RetentionTime ?? 0.0;
                return aKey.CompareTo(bKey);
            });
            int n = sorted.Count;

            var decilesSeen = new HashSet<int>();
            for (int i = 0; i < n; i++)
            {
                if (!selSet.Contains(sorted[i].Sequence)) continue;
                int d = (int)((long)i * deciles / n);
                d = Math.Min(d, deciles - 1);
                decilesSeen.Add(d);
            }
            return decilesSeen.Count;
        }

        // ────────────────────────────────────────────────────────────────
        //  Key helpers
        // ────────────────────────────────────────────────────────────────

        // Key that uniquely identifies a benchmark precursor (sequence + charge)
        private static string PrecursorKey(LibraryPrecursorInput p)
            => $"{p.Sequence}/{p.ChargeState}";

        // Key built from a DiaSearchResult — same format
        private static string ResultKey(DiaSearchResult r)
            => $"{r.Sequence}/{r.ChargeState}";

        // ────────────────────────────────────────────────────────────────
        //  Statistics helpers
        // ────────────────────────────────────────────────────────────────

        private static float ComputeMedian(List<float> values)
        {
            if (values.Count == 0) return float.NaN;
            var sorted = new List<float>(values);
            sorted.Sort();
            int mid = sorted.Count / 2;
            return sorted.Count % 2 == 0
                ? (sorted[mid - 1] + sorted[mid]) * 0.5f
                : sorted[mid];
        }

        // ────────────────────────────────────────────────────────────────
        //  Console output helpers
        // ────────────────────────────────────────────────────────────────

        private static void PrintAccuracyTable(WindowSummary[] summaries)
        {
            Console.WriteLine("╔══════════════════════════════════════════════════════════════════════╗");
            Console.WriteLine("║          Peak Group Selection Accuracy vs Window Width               ║");
            Console.WriteLine("╠══════════╦══════════════╦═══════════════╦═════════════╦═════════════╣");
            Console.WriteLine("║ Window   ║ Correct      ║ Marginal      ║ Wrong       ║ MedianErr   ║");
            Console.WriteLine("║          ║ (<0.5 min)   ║ (0.5–1.5 min) ║ (>1.5 min)  ║             ║");
            Console.WriteLine("╠══════════╬══════════════╬═══════════════╬═════════════╬═════════════╣");
            foreach (var s in summaries)
            {
                int n = s.Total;
                string win = $"±{s.HalfWidth:F0} min".PadRight(8);
                string cor = FormatCell(s.Correct, n, 12);
                string marg = FormatCell(s.Marginal, n, 13);
                string wrg = FormatCell(s.Wrong, n, 11);
                string merr = float.IsNaN(s.MedianError)
                    ? "  N/A      "
                    : $"  {s.MedianError:F3} min".PadRight(13);
                Console.WriteLine($"║ {win} ║{cor}║{marg}║{wrg}║{merr}║");
            }
            Console.WriteLine("╚══════════╩══════════════╩═══════════════╩═════════════╩═════════════╝");
            Console.WriteLine();
        }

        private static string FormatCell(int count, int total, int width)
        {
            double pct = total > 0 ? 100.0 * count / total : 0.0;
            return $"  {count}/{total} {pct:F0}%".PadRight(width);
        }

        private static void PrintWorstMisses(List<PeakAccuracyRecord> allRecords, float worstWindowHalf, int topN)
        {
            var windowRecords = new List<PeakAccuracyRecord>();
            for (int i = 0; i < allRecords.Count; i++)
                if (Math.Abs(allRecords[i].WindowHalf - worstWindowHalf) < 0.01f)
                    windowRecords.Add(allRecords[i]);

            if (windowRecords.Count == 0)
            {
                Console.WriteLine($"No records for ±{worstWindowHalf:F0} min window.");
                return;
            }

            windowRecords.Sort((a, b) => b.AbsError.CompareTo(a.AbsError));

            Console.WriteLine($"Top {topN} worst peak selection misses at ±{worstWindowHalf:F0} min window:");
            Console.WriteLine("  Rank | Sequence            | iRT   | TrueRT | SelectedApex | Error  | Cosine");
            Console.WriteLine("  ─────┼─────────────────────┼───────┼────────┼──────────────┼────────┼───────");

            int shown = 0;
            for (int i = 0; i < windowRecords.Count && shown < topN; i++)
            {
                var r = windowRecords[i];
                // Skip "correct" entries — we want actual misses only
                if (r.Classification == Correct) continue;

                string seqStr = $"{TruncateSeq(r.Sequence)}/{r.ChargeState}".PadRight(20);
                string irt = float.IsNaN(r.IRt) ? "  N/A " : r.IRt.ToString("F1").PadLeft(6);
                string trueRt = r.TrueRt.ToString("F1").PadLeft(6);
                string selApex = float.IsNaN(r.SelectedApexRt) ? "      NaN" : r.SelectedApexRt.ToString("F1").PadLeft(9);
                string err = (r.AbsError >= 99f ? ">5.0" : r.AbsError.ToString("F2")).PadLeft(6);
                string cos = float.IsNaN(r.ApexCosine) ? "  N/A" : r.ApexCosine.ToString("F3").PadLeft(5);

                Console.WriteLine($"  {shown + 1,4} | {seqStr}| {irt} | {trueRt} | {selApex}   | {err}   | {cos}");
                shown++;
            }
            Console.WriteLine();
        }

        private static string TruncateSeq(string seq)
        {
            return seq.Length > 18 ? seq.Substring(0, 17) + "…" : seq;
        }

        // ────────────────────────────────────────────────────────────────
        //  TSV export
        // ────────────────────────────────────────────────────────────────

        private static void WriteDetailTsv(List<PeakAccuracyRecord> records, string path)
        {
            // Ensure output directory exists
            string dir = Path.GetDirectoryName(path);
            if (!string.IsNullOrEmpty(dir) && !Directory.Exists(dir))
                Directory.CreateDirectory(dir);

            using var w = new StreamWriter(path);
            w.WriteLine("Sequence\tChargeState\tiRT\tTrueRT\tWindowHalf\tSelectedApexRT\tAbsError\tClassification\tApexCosine\tCandidateCount");

            for (int i = 0; i < records.Count; i++)
            {
                var r = records[i];
                w.Write(r.Sequence);
                w.Write('\t'); w.Write(r.ChargeState);
                w.Write('\t'); w.Write(float.IsNaN(r.IRt) ? "NA" : r.IRt.ToString("F4"));
                w.Write('\t'); w.Write(r.TrueRt.ToString("F4"));
                w.Write('\t'); w.Write(r.WindowHalf.ToString("F1"));
                w.Write('\t'); w.Write(float.IsNaN(r.SelectedApexRt) ? "NaN" : r.SelectedApexRt.ToString("F4"));
                w.Write('\t'); w.Write(r.AbsError.ToString("F4"));
                w.Write('\t'); w.Write(r.Classification);
                w.Write('\t'); w.Write(float.IsNaN(r.ApexCosine) ? "NaN" : r.ApexCosine.ToString("F4"));
                w.Write('\t'); w.Write(r.CandidateCount);
                w.WriteLine();
            }
        }

        // ────────────────────────────────────────────────────────────────
        //  Internal data types
        // ────────────────────────────────────────────────────────────────

        private struct PeakAccuracyRecord
        {
            public string Sequence;
            public int ChargeState;
            public float IRt;
            public float TrueRt;
            public float WindowHalf;
            public float SelectedApexRt;
            public float AbsError;
            public string Classification;
            public float ApexCosine;
            public int CandidateCount;
        }

        private struct WindowSummary
        {
            public float HalfWidth;
            public int Total;
            public int Correct;
            public int Marginal;
            public int Wrong;
            public float MedianError;
        }
    }
}