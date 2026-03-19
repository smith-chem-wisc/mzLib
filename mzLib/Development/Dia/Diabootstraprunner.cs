// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// Location: mzLib/Development/Dia/DiaBootstrapRunner.cs
//
// Development runner for the bidirectional bootstrap calibrator.
// All calibration logic lives in DiaBootstrapCalibrator.RunBidirectional().
// This file contains only file loading and display scaffolding.

using Easy.Common.Extensions;
using MassSpectrometry;
using MassSpectrometry.Dia;
using MassSpectrometry.Dia.Calibration;
using Readers;
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace Development.Dia
{
    public static class DiaBootstrapRunner
    {
        public static void Run(
            string rawFilePath,
            string[] mspPaths,
            float ppmTol = 20f,
            bool filenameDecoySuffix = true)
        {
            Header("DIA Bootstrap Calibrator — Development Runner");

            Section("1. Loading raw file");
            MsDataScan[]? allScans = LoadRawFile(rawFilePath);
            if (allScans == null) return;

            Section("2. Building DIA scan index");
            DiaScanIndex? scanIndex = BuildIndex(allScans);
            if (scanIndex == null) return;

            Section("3. Loading MSP libraries");
            List<LibraryPrecursorInput>? precursors = LoadMspLibraries(mspPaths, filenameDecoySuffix);
            if (precursors == null || precursors.Count == 0)
            { Error("No library precursors loaded."); return; }
            PrintLibrarySummary(precursors);

            Section("4. Bidirectional sliding-window bootstrap");
            var totalSw = Stopwatch.StartNew();

            BootstrapCalibrationResult result = DiaBootstrapCalibrator.RunBidirectional(
                precursors, scanIndex,
                ppmTolerance: ppmTol,
                progressReporter: msg => Console.WriteLine("  " + msg));

            totalSw.Stop();
            Console.WriteLine($"\n  Total bootstrap time: {totalSw.Elapsed.TotalSeconds:F1} s");

            if (result.Seed == null)
            { Error("Bootstrap failed at seed step."); return; }

            Section("5. Per-step summary");
            PrintStepSummary(result, scanIndex);

            Section("6. Global calibration");
            if (result.GlobalOlsLine != null)
            {
                double cov2 = Coverage(result.GlobalOlsLine.ResidualSigma, 2f);
                Console.WriteLine($"  OLS line:  {result.GlobalOlsLine}");
                Console.WriteLine($"  OLS 2σ hw: {result.GlobalOlsLine.ResidualSigma * 2f:F4} min" +
                                  $"  → {cov2 * 100:F1}% coverage");
            }
            else
            {
                Warn("Global OLS fit failed.");
            }

            Console.WriteLine();

            if (result.LowessModel != null)
            {
                Console.WriteLine($"  LOWESS: R²={result.LowessModel.RSquared:F4}" +
                                  $"  σ={result.LowessModel.SigmaMinutes:F4} min" +
                                  $"  n={result.LowessModel.AnchorCount}");
                Console.WriteLine($"  LOWESS 2σ hw: {result.LowessModel.SigmaMinutes * 2:F4} min");
                Console.WriteLine();
                PrintLowessDiagnostics(result.LowessModel);
            }
            else
            {
                Warn("LOWESS fit failed (insufficient anchors).");
            }
        }

        // ══════════════════════════════════════════════════════════════════
        // Display
        // ══════════════════════════════════════════════════════════════════

        private static void PrintStepSummary(
            BootstrapCalibrationResult result, DiaScanIndex scanIndex)
        {
            Console.WriteLine();
            Console.WriteLine($"  {"Arm",-6}  {"Range",-28}  {"Sep",5}  {"R²",6}  {"σ min",6}  {"hw min",6}  {"Anchors",7}");
            Console.WriteLine("  " + new string('─', 76));

            var all = new List<(string Tag, BootstrapSliceResult R)>();
            all.Add(("SEED ", result.Seed!));
            foreach (var s in result.LeftSteps) all.Add(("LEFT ", s));
            foreach (var s in result.RightSteps) all.Add(("RIGHT", s));
            all = all.OrderBy(x => x.R.SliceLo).ToList();

            foreach (var (tag, r) in all)
            {
                string sep = $"{r.Model.SeparationRatio:F2}";
                string r2 = r.RtLine != null ? $"{r.RtLine.RSquared:F4}" : "  n/a";
                string sig = r.RtLine != null ? $"{r.RtLine.ResidualSigma:F4}" : "  n/a";
                string hw = r.RtLine != null ? $"{r.RtLine.HalfWidth:F4}" : "  n/a";
                string pct = ToPercentileLabel(scanIndex, r.SliceLo, r.SliceHi);
                Console.WriteLine(
                    $"  {tag}  {pct,-28}  {sep,5}  {r2,6}  {sig,6}  {hw,6}  {r.Anchors.Count,7}");
            }

            Console.WriteLine();
            Console.WriteLine($"  Total anchors: {result.AllAnchors.Count}" +
                $"  (seed={result.Seed!.Anchors.Count}" +
                $"  left={result.LeftSteps.Sum(s => s.Anchors.Count)}" +
                $"  right={result.RightSteps.Sum(s => s.Anchors.Count)})");
        }

        private static void PrintLowessDiagnostics(LowessRtModel m)
        {
            Console.WriteLine("  LOWESS curve (library RT → predicted obs RT, with local σ):");
            Console.WriteLine($"  {"LibRT",8}  {"PredRT",8}  {"LocalSigma",12}  {"hw 2σ",8}");
            Console.WriteLine("  " + new string('-', 46));

            int n = m.FittedLibraryRts.Length;
            int nSamples = 20;
            for (int s = 0; s < nSamples; s++)
            {
                int idx = (int)((double)s / (nSamples - 1) * (n - 1));
                double lib = m.FittedLibraryRts[idx];
                double obs = m.ToMinutes(lib);
                double sig = m.GetLocalSigma(lib);
                Console.WriteLine($"  {lib,8:F2}  {obs,8:F3}  {sig,12:F4}  {sig * 2,8:F4}");
            }

            Console.WriteLine();
            Console.WriteLine($"  Global LOWESS σ:  {m.SigmaMinutes:F4} min");
            Console.WriteLine($"  Global LOWESS R²: {m.RSquared:F4}");
            Console.WriteLine($"  Monotonic enforced: {m.MonotonicEnforced}");
            Console.WriteLine();
            Console.WriteLine("  Usage in full-run search:");
            Console.WriteLine("    predictedRt = model.ToMinutes(libraryRt)");
            Console.WriteLine("    hw = Math.Min(model.GetLocalSigma(libraryRt), 1.5) * 2.0");
        }

        // ══════════════════════════════════════════════════════════════════
        // File loading
        // ══════════════════════════════════════════════════════════════════

        private static MsDataScan[]? LoadRawFile(string path)
        {
            Console.WriteLine($"  Path:   {path}");
            Console.WriteLine($"  Format: {Path.GetExtension(path).ToUpperInvariant()}");
            if (!File.Exists(path)) { Error($"File not found: {path}"); return null; }
            var sw = Stopwatch.StartNew();
            try
            {
                var file = MsDataFileReader.GetDataFile(path);
                file.LoadAllStaticData();
                var scans = file.GetAllScansList().ToArray();
                sw.Stop();
                Console.WriteLine($"  Loaded {scans.Length} scans in {sw.Elapsed.TotalSeconds:F2} s" +
                    $"  (MS1={scans.Count(s => s.MsnOrder == 1)}," +
                    $" MS2={scans.Count(s => s.MsnOrder == 2)}," +
                    $" RT {scans.First().RetentionTime:F2}–{scans.Last().RetentionTime:F2} min)");
                return scans;
            }
            catch (Exception ex) { Error($"Load failed: {ex.Message}"); return null; }
        }

        private static DiaScanIndex? BuildIndex(MsDataScan[] scans)
        {
            var sw = Stopwatch.StartNew();
            try
            {
                var idx = DiaScanIndexBuilder.Build(scans);
                sw.Stop();
                (float lo, float hi) = DiaBootstrapCalibrator.ComputeSliceBounds(idx);
                Console.WriteLine($"  Built in {sw.ElapsedMilliseconds:F0} ms" +
                    $"  ({idx.ScanCount:N0} MS2, {idx.TotalPeakCount:N0} peaks," +
                    $" {idx.WindowCount} windows)");
                Console.WriteLine($"  Seed slice: [{lo:F3}, {hi:F3}] min");
                return idx;
            }
            catch (Exception ex) { Error($"Failed: {ex.Message}"); return null; }
        }

        private static List<LibraryPrecursorInput>? LoadMspLibraries(
            string[] paths, bool filenameDecoySuffix)
        {
            var all = new List<LibraryPrecursorInput>();
            foreach (string p in paths)
            {
                if (!File.Exists(p)) { Warn($"MSP not found: {p}"); continue; }
                bool forceDecoy = filenameDecoySuffix &&
                    Path.GetFileNameWithoutExtension(p)
                        .IndexOf("decoy", StringComparison.OrdinalIgnoreCase) >= 0;
                Console.WriteLine($"  {Path.GetFileName(p)}  [{(forceDecoy ? "DECOY" : "TARGET")}]");
                var sw = Stopwatch.StartNew();
                var entries = KoinaMspParser.Parse(p);
                sw.Stop();
                if (forceDecoy)
                {
                    var ov = new List<LibraryPrecursorInput>(entries.Count);
                    foreach (var e in entries)
                        ov.Add(new LibraryPrecursorInput(e.Sequence, e.PrecursorMz, e.ChargeState,
                            e.RetentionTime, true, e.FragmentMzs, e.FragmentIntensities, e.IrtValue));
                    entries = ov;
                }
                int t = entries.Count(e => !e.IsDecoy), d = entries.Count(e => e.IsDecoy);
                Console.WriteLine($"    {entries.Count} entries ({t} T, {d} D)" +
                                  $" in {sw.ElapsedMilliseconds:F0} ms");
                all.AddRange(entries);
            }
            return all;
        }

        private static void PrintLibrarySummary(IReadOnlyList<LibraryPrecursorInput> prec)
        {
            int t = prec.Count(p => !p.IsDecoy), d = prec.Count(p => p.IsDecoy);
            int hasRt = prec.Count(p => p.RetentionTime.HasValue);
            Console.WriteLine($"  Total: {prec.Count} ({t} T, {d} D)" +
                $"  avg frags={prec.Average(p => p.FragmentCount):F1}" +
                $"  with RT={hasRt} ({hasRt * 100.0 / prec.Count:F1}%)");
        }

        // ══════════════════════════════════════════════════════════════════
        // Formatting utilities
        // ══════════════════════════════════════════════════════════════════

        private static string ToPercentileLabel(DiaScanIndex idx, float lo, float hi)
        {
            int n = idx.ScanCount;
            if (n == 0) return $"[{lo:F1},{hi:F1}]";
            float[] buf = ArrayPool<float>.Shared.Rent(n);
            try
            {
                for (int i = 0; i < n; i++) buf[i] = idx.GetScanRt(i);
                Array.Sort(buf, 0, n);
                float range = buf[n - 1] - buf[0];
                if (range <= 0) return $"[{lo:F1},{hi:F1}]";
                int loPct = Math.Clamp((int)Math.Round((lo - buf[0]) / range * 100), 0, 100);
                int hiPct = Math.Clamp((int)Math.Round((hi - buf[0]) / range * 100), 0, 100);
                return $"{loPct}–{hiPct}%  [{lo:F2},{hi:F2}]min";
            }
            finally { ArrayPool<float>.Shared.Return(buf); }
        }

        private static double Coverage(float sigma, float sigmaMultiplier)
        {
            double z = sigmaMultiplier;
            double t = 1.0 / (1.0 + 0.3275911 * z);
            double poly = t * (0.254829592 + t * (-0.284496736
                        + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
            return 1.0 - poly * Math.Exp(-z * z);
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

        private static void Error(string m) => Console.WriteLine("  [ERROR] " + m);
        private static void Warn(string m) => Console.WriteLine("  [WARN]  " + m);
    }
}