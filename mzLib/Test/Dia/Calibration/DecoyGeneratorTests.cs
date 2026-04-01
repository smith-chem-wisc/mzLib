// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// Location: Test/Dia/DecoyGeneratorTests.cs
//
// NUnit 4 test that:
//   1. Loads the DIA-NN MSP library
//   2. Generates decoys via DecoyGenerator
//   3. Compares each paired target/decoy
//   4. Writes a detailed diagnostic TSV to F:\DiaBenchmark\NewLibrary\

using MassSpectrometry.Dia;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Test.Dia
{
    [TestFixture]
    public class DecoyGeneratorTests
    {
        // ── Paths ─────────────────────────────────────────────────────────────
        private const string MspPath =
            @"F:\DiaBenchmark\PXD005573\DiannOut\diannFig2helalib.msp";

        private const string OutputDir =
            @"F:\DiaBenchmark\NewLibrary";

        private const string PairTsvPath =
            @"F:\DiaBenchmark\NewLibrary\decoy_pair_diagnostics.tsv";

        private const string SummaryPath =
            @"F:\DiaBenchmark\NewLibrary\decoy_summary.txt";

        // ── Test ──────────────────────────────────────────────────────────────

        [Test]
        public void DecoyPairDiagnostics()
        {
            // ── 1. Setup ──────────────────────────────────────────────────────
            Assert.That(File.Exists(MspPath),
                $"MSP library not found: {MspPath}");

            Directory.CreateDirectory(OutputDir);

            // ── 2. Load targets ───────────────────────────────────────────────
            Console.WriteLine($"Loading library: {MspPath}");
            var targets = KoinaMspParser.Parse(MspPath, minIntensity: 0.05f);
            Console.WriteLine($"Loaded {targets.Count:N0} targets");
            Assert.That(targets.Count, Is.GreaterThan(0), "No targets loaded");

            // ── 3. Generate decoys ────────────────────────────────────────────
            Console.WriteLine("Generating decoys...");
            var decoys = DecoyGenerator.GenerateFromTargets(targets, randomSeed: 42);
            Assert.That(decoys.Count, Is.EqualTo(targets.Count),
                "Decoy count must equal target count (1:1 ratio)");

            // ── 4. Per-pair analysis ──────────────────────────────────────────
            var inv = CultureInfo.InvariantCulture;

            // Counters for summary
            int palindromes = 0;  // scrambled == forward
            int shortPeptides = 0;  // length <= 2
            int zeroIntDecoys = 0;  // decoy has any zero-intensity fragment
            int mzOverlapCount = 0;  // any decoy frag within 20 ppm of target frag
            int totalFragOverlaps = 0;  // total fragment m/z overlaps across all pairs
            int totalTargetFrags = 0;
            int totalDecoyFrags = 0;
            var intensityCorrs = new List<double>();  // Pearson r of intensities
            var mzShifts = new List<double>();  // mean |decoy_mz - nearest_target_mz|

            // Intensity overlap metric: how similar are the intensity vectors?
            // We compare decoy intensities vs target intensities at same rank
            // (both sorted descending by intensity).
            var intensityOverlapRatios = new List<double>();

            Console.WriteLine($"Writing pair diagnostics to: {PairTsvPath}");
            var w = new StreamWriter(PairTsvPath);

            // Header
            w.WriteLine(string.Join("\t", new[]
            {
                "Index",
                "TargetSequence",
                "DecoySequence",
                "Charge",
                "PrecursorMz",
                "PeptideLength",
                "TargetFragCount",
                "DecoyFragCount",
                "IsPalindrome",
                "MzOverlapCount",          // decoy frags within 20ppm of any target frag
                "MzOverlapFraction",       // fraction of decoy frags that overlap
                "MeanDecoyMzShift",        // mean |decoy_mz - nearest_target_mz| in Da
                "MaxDecoyMzShift",
                "NDP",                     // normalized dot product in m/z space (20 ppm tolerance)
                "TopIntensityRatio",       // decoy top intensity / target top intensity
                "HasZeroIntDecoyFrag",
                "TargetRt",
                "DecoyRt",
                "TargetTopMz",
                "DecoyTopMz",
                "TargetTopInt",
                "DecoyTopInt",
                "TargetMzRange",           // max-min of target fragment mzs
                "DecoyMzRange",
            }));

            for (int i = 0; i < targets.Count; i++)
            {
                var t = targets[i];
                var d = decoys[i];

                string tSeq = StripMods(t.Sequence ?? "");
                string dSeq = StripMods(d.Sequence ?? "");
                int pepLen = tSeq.Length;

                bool isPalin = tSeq == dSeq;
                if (isPalin) palindromes++;
                if (pepLen <= 2) shortPeptides++;

                int tFragCount = t.FragmentMzs?.Length ?? 0;
                int dFragCount = d.FragmentMzs?.Length ?? 0;
                totalTargetFrags += tFragCount;
                totalDecoyFrags += dFragCount;

                // ── m/z overlap (within 20 ppm) ──────────────────────────────
                int overlapCount = 0;
                var mzShiftList = new List<double>(dFragCount);
                bool hasZeroInt = false;

                if (tFragCount > 0 && dFragCount > 0)
                {
                    for (int j = 0; j < dFragCount; j++)
                    {
                        float dMz = d.FragmentMzs[j];
                        float dInt = d.FragmentIntensities?[j] ?? 0f;
                        if (dInt == 0f) hasZeroInt = true;

                        // Find nearest target fragment
                        float nearest = FindNearest(t.FragmentMzs, dMz);
                        double shift = Math.Abs(dMz - nearest);
                        mzShiftList.Add(shift);

                        double ppm = shift / nearest * 1e6;
                        if (ppm <= 20.0) overlapCount++;
                    }
                }

                if (hasZeroInt) zeroIntDecoys++;
                if (overlapCount > 0) mzOverlapCount++;
                totalFragOverlaps += overlapCount;

                double meanShift = mzShiftList.Count > 0 ? mzShiftList.Average() : double.NaN;
                double maxShift = mzShiftList.Count > 0 ? mzShiftList.Max() : double.NaN;
                double overlapFrac = dFragCount > 0 ? (double)overlapCount / dFragCount : 0.0;

                mzShifts.Add(double.IsNaN(meanShift) ? 0 : meanShift);

                // ── Spectral similarity: normalized dot product ───────────────
                // Build sparse intensity vectors aligned to a common m/z grid.
                // Two fragments are considered the same peak if within 20 ppm.
                // NDP = dot(T,D) / (|T| * |D|)  where vectors are in m/z space.
                // NDP=1 means identical spectra, NDP=0 means no shared peaks.
                double intCorr = double.NaN;   // reused as NDP for output column
                double topIntRatio = double.NaN;

                if (tFragCount > 0 && dFragCount > 0 &&
                    t.FragmentMzs != null && d.FragmentMzs != null &&
                    t.FragmentIntensities != null && d.FragmentIntensities != null)
                {
                    intCorr = NormalizedDotProduct(
                        t.FragmentMzs, t.FragmentIntensities,
                        d.FragmentMzs, d.FragmentIntensities,
                        ppmTolerance: 20.0);

                    float tMaxInt = t.FragmentIntensities.Max();
                    float dMaxInt = d.FragmentIntensities.Max();
                    topIntRatio = tMaxInt > 0 ? dMaxInt / tMaxInt : double.NaN;
                }

                if (!double.IsNaN(intCorr)) intensityCorrs.Add(intCorr);
                if (!double.IsNaN(topIntRatio)) intensityOverlapRatios.Add(topIntRatio);

                // Target/decoy top fragment info
                float tTopMz = tFragCount > 0 ? t.FragmentMzs.Max() : float.NaN;
                float dTopMz = dFragCount > 0 ? d.FragmentMzs.Max() : float.NaN;
                float tTopInt = (tFragCount > 0 && t.FragmentIntensities != null)
                    ? t.FragmentIntensities.Max() : float.NaN;
                float dTopInt = (dFragCount > 0 && d.FragmentIntensities != null)
                    ? d.FragmentIntensities.Max() : float.NaN;
                float tMzRange = tFragCount > 1
                    ? t.FragmentMzs.Max() - t.FragmentMzs.Min() : float.NaN;
                float dMzRange = dFragCount > 1
                    ? d.FragmentMzs.Max() - d.FragmentMzs.Min() : float.NaN;

                // ── Write row ─────────────────────────────────────────────────
                w.WriteLine(string.Join("\t", new[]
                {
                    i.ToString(inv),
                    t.Sequence ?? "",
                    d.Sequence ?? "",
                    t.ChargeState.ToString(inv),
                    t.PrecursorMz.ToString("F5", inv),
                    pepLen.ToString(inv),
                    tFragCount.ToString(inv),
                    dFragCount.ToString(inv),
                    isPalin ? "1" : "0",
                    overlapCount.ToString(inv),
                    overlapFrac.ToString("F4", inv),
                    double.IsNaN(meanShift) ? "NA" : meanShift.ToString("F4", inv),
                    double.IsNaN(maxShift)  ? "NA" : maxShift.ToString("F4", inv),
                    double.IsNaN(intCorr)   ? "NA" : intCorr.ToString("F4", inv),
                    double.IsNaN(topIntRatio) ? "NA" : topIntRatio.ToString("F4", inv),
                    hasZeroInt ? "1" : "0",
                    t.RetentionTime.HasValue ? t.RetentionTime.Value.ToString("F4", inv) : "NA",
                    d.RetentionTime.HasValue ? d.RetentionTime.Value.ToString("F4", inv) : "NA",
                    float.IsNaN(tTopMz) ? "NA" : tTopMz.ToString("F4", inv),
                    float.IsNaN(dTopMz) ? "NA" : dTopMz.ToString("F4", inv),
                    float.IsNaN(tTopInt) ? "NA" : tTopInt.ToString("F4", inv),
                    float.IsNaN(dTopInt) ? "NA" : dTopInt.ToString("F4", inv),
                    float.IsNaN(tMzRange) ? "NA" : tMzRange.ToString("F4", inv),
                    float.IsNaN(dMzRange) ? "NA" : dMzRange.ToString("F4", inv),
                }));
            }

            // Close writer before reading back for CountHighOverlap
            w.Close();
            w.Dispose();

            // ── 4b. RT diagnostic ─────────────────────────────────────────────
            int targetNoRt = 0, decoyNoRt = 0;
            int targetHasRt = 0, decoyHasRt = 0;
            for (int i = 0; i < targets.Count; i++)
            {
                if (targets[i].RetentionTime.HasValue) targetHasRt++; else targetNoRt++;
                if (decoys[i].RetentionTime.HasValue) decoyHasRt++; else decoyNoRt++;
            }

            // ── 5. Summary ────────────────────────────────────────────────────
            double meanCorr = intensityCorrs.Count > 0 ? intensityCorrs.Average() : double.NaN;
            double meanMzShift = mzShifts.Count > 0 ? mzShifts.Average() : double.NaN;
            double meanTopRatio = intensityOverlapRatios.Count > 0
                ? intensityOverlapRatios.Average() : double.NaN;

            // Fraction of decoy frags that overlap with target frags within 20 ppm
            double overallOverlapFrac = totalDecoyFrags > 0
                ? (double)totalFragOverlaps / totalDecoyFrags : 0.0;

            var summary = new List<string>
            {
                "═══════════════════════════════════════════════════════════",
                " Decoy Generator Diagnostics",
                "═══════════════════════════════════════════════════════════",
                $" Library:              {Path.GetFileName(MspPath)}",
                $" Target precursors:    {targets.Count:N0}",
                $" Decoy precursors:     {decoys.Count:N0}",
                "",
                "── Sequence Quality ───────────────────────────────────────",
                $" Palindromes (decoy == target):  {palindromes:N0} ({palindromes * 100.0 / targets.Count:F2}%)",
                $" Short peptides (len ≤ 2):       {shortPeptides:N0}",
                "",
                "── Fragment m/z Overlap (within 20 ppm) ──────────────────",
                $" Pairs with ≥1 overlapping frag: {mzOverlapCount:N0} ({mzOverlapCount * 100.0 / targets.Count:F2}%)",
                $" Total overlapping frags:        {totalFragOverlaps:N0} / {totalDecoyFrags:N0} decoy frags",
                $" Overall overlap fraction:       {overallOverlapFrac:F4} ({overallOverlapFrac * 100:F2}%)",
                $" Mean |decoy_mz - nearest_target_mz|: {meanMzShift:F4} Da",
                "",
                "── Intensity Patterns ─────────────────────────────────────",
                $" Mean spectral similarity (NDP, 20 ppm): {meanCorr:F4}",
                $"   (0=no shared peaks, 1=identical spectra)",
                $" Mean top-intensity ratio (decoy/target):     {meanTopRatio:F4}",
                $"   (should be ≈1.0 since both come from same source)",
                $" Pairs with any zero-intensity decoy frag:   {zeroIntDecoys:N0}",
                "",
                "── Retention Time Coverage ───────────────────────────────",
                $" Targets with RT:    {targetHasRt:N0} / {targets.Count:N0} ({targetHasRt * 100.0 / targets.Count:F2}%)",
                $" Targets without RT: {targetNoRt:N0} ({targetNoRt * 100.0 / targets.Count:F2}%)",
                $" Decoys with RT:     {decoyHasRt:N0} / {decoys.Count:N0} ({decoyHasRt * 100.0 / decoys.Count:F2}%)",
                $" Decoys without RT:  {decoyNoRt:N0} ({decoyNoRt * 100.0 / decoys.Count:F2}%)",
                "",
                "── Fragment Counts ────────────────────────────────────────",
                $" Mean target frags/precursor: {(targets.Count > 0 ? (double)totalTargetFrags / targets.Count : 0):F2}",
                $" Mean decoy frags/precursor:  {(decoys.Count > 0 ? (double)totalDecoyFrags / decoys.Count : 0):F2}",
                "",
                "── Key Concern Indicators ─────────────────────────────────",
                $" HIGH spectral similarity NDP >0.5: {intensityCorrs.Count(c => c > 0.5):N0} pairs ({intensityCorrs.Count(c => c > 0.5) * 100.0 / Math.Max(intensityCorrs.Count, 1):F1}%)",
                $"   → These decoys share significant signal with their targets (bad null dist)",
                $" HIGH spectral similarity NDP >0.8: {intensityCorrs.Count(c => c > 0.8):N0} pairs ({intensityCorrs.Count(c => c > 0.8) * 100.0 / Math.Max(intensityCorrs.Count, 1):F1}%)",
                $"   → These decoys are nearly identical to their targets (very bad)",
                $" HIGH m/z overlap (>50% frags):    {CountHighOverlap(PairTsvPath):N0} pairs",
                $"   → These decoys share fragment m/z with targets (bad null dist)",
                "",
                "── Output Files ───────────────────────────────────────────",
                $" Pair diagnostics TSV: {PairTsvPath}",
                $" This summary:         {SummaryPath}",
                "═══════════════════════════════════════════════════════════",
            };

            foreach (var line in summary) Console.WriteLine(line);
            File.WriteAllLines(SummaryPath, summary);

            // ── 6. Assertions ─────────────────────────────────────────────────
            Assert.That(decoys.Count, Is.EqualTo(targets.Count),
                "1:1 target:decoy ratio required");

            Assert.That(targetNoRt, Is.EqualTo(0),
                $"{targetNoRt} targets have no RetentionTime — every library spectrum must have an RT");

            Assert.That(decoyNoRt, Is.EqualTo(0),
                $"{decoyNoRt} decoys have no RetentionTime — decoys inherit RT from their target");

            Assert.That(palindromes * 100.0 / targets.Count, Is.LessThan(5.0),
                $"Too many palindromes: {palindromes} ({palindromes * 100.0 / targets.Count:F1}%)");

            Assert.That(overallOverlapFrac, Is.LessThan(0.30),
                $"Too many fragment m/z overlaps: {overallOverlapFrac:F3} (threshold 0.30) — " +
                "decoys share too many fragment m/z values with targets within 20 ppm");

            // Warn (not fail) if spectral similarity is high
            if (!double.IsNaN(meanCorr) && meanCorr > 0.1)
                Console.WriteLine($"  WARNING: Mean NDP = {meanCorr:F4} > 0.1 " +
                    "— decoys share more signal with targets than expected for random sequences");

            Console.WriteLine($"\nDiagnostic files written to: {OutputDir}");
        }

        // ── Helpers ───────────────────────────────────────────────────────────

        private static string StripMods(string seq)
        {
            var sb = new System.Text.StringBuilder(seq.Length);
            bool inMod = false;
            char open = ' ';
            foreach (char c in seq)
            {
                if (!inMod)
                {
                    if (c == '(' || c == '[') { inMod = true; open = c; }
                    else sb.Append(c);
                }
                else if (c == (open == '(' ? ')' : ']'))
                {
                    inMod = false;
                }
            }
            return sb.ToString();
        }

        private static float FindNearest(float[] arr, float target)
        {
            if (arr == null || arr.Length == 0) return float.NaN;
            int lo = 0, hi = arr.Length - 1;
            while (lo < hi)
            {
                int mid = (lo + hi) / 2;
                if (arr[mid] < target) lo = mid + 1; else hi = mid;
            }
            if (lo == 0) return arr[0];
            return Math.Abs(arr[lo] - target) < Math.Abs(arr[lo - 1] - target)
                ? arr[lo] : arr[lo - 1];
        }

        /// <summary>
        /// Normalized dot product between two spectra in m/z space.
        /// Fragments are matched if within ppmTolerance ppm of each other.
        /// NDP = dot(T,D) / (norm(T) * norm(D))  range [0,1].
        /// NDP=0: no shared peaks. NDP=1: identical spectra.
        /// </summary>
        private static double NormalizedDotProduct(
            float[] tMzs, float[] tInts,
            float[] dMzs, float[] dInts,
            double ppmTolerance = 20.0)
        {
            double dot = 0.0, normT = 0.0, normD = 0.0;

            // Accumulate norms
            for (int i = 0; i < tInts.Length; i++) normT += tInts[i] * tInts[i];
            for (int j = 0; j < dInts.Length; j++) normD += dInts[j] * dInts[j];

            normT = Math.Sqrt(normT);
            normD = Math.Sqrt(normD);
            if (normT < 1e-12 || normD < 1e-12) return 0.0;

            // Two-pointer walk over sorted m/z arrays to find matching peaks
            int ti = 0, di = 0;
            while (ti < tMzs.Length && di < dMzs.Length)
            {
                double ppm = Math.Abs(tMzs[ti] - dMzs[di]) / tMzs[ti] * 1e6;
                if (ppm <= ppmTolerance)
                {
                    dot += tInts[ti] * dInts[di];
                    ti++; di++;
                }
                else if (tMzs[ti] < dMzs[di]) ti++;
                else di++;
            }

            return dot / (normT * normD);
        }

        /// <summary>
        /// Counts pairs from the already-written TSV where MzOverlapFraction > 0.5.
        /// Reads back the file to avoid recomputing.
        /// </summary>
        private static int CountHighOverlap(string tsvPath)
        {
            if (!File.Exists(tsvPath)) return 0;
            int col = -1, count = 0;
            bool header = true;
            foreach (var line in File.ReadLines(tsvPath))
            {
                var f = line.Split('\t');
                if (header)
                {
                    for (int i = 0; i < f.Length; i++)
                        if (f[i] == "MzOverlapFraction") { col = i; break; }
                    header = false; continue;
                }
                if (col >= 0 && col < f.Length &&
                    double.TryParse(f[col], NumberStyles.Float,
                        CultureInfo.InvariantCulture, out double v) && v > 0.5)
                    count++;
            }
            return count;
        }
    }
}