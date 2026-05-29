using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Readers;

namespace Test.MassSpectrometryTests.Deconvolution
{
    /// <summary>
    /// EXPLICIT per-spectrum fidelity test: for the SAME real MS1 scans, compare MetaFlashDecon's
    /// deconvolved monoisotopic masses against real OpenMS FLASHDeconv's (parsed from the per-scan
    /// _ms1.msalign ground truth). Per scan it reports matched / ±1-Da-isotope-off / our-extra /
    /// real-missed masses — directly identifying where our per-spectrum decon diverges from real
    /// FLASHDeconv on the actual spectra we get wrong (the user's idea).
    ///
    ///   dotnet test --filter FullyQualifiedName~MetaFlashDeconVsRealPerSpectrum
    /// </summary>
    [TestFixture]
    [Explicit("Per-spectrum MetaFlashDecon vs real FLASHDeconv (_ms1.msalign); needs E:\\Projects\\JurkatTopDown.")]
    [ExcludeFromCodeCoverage]
    public class MetaFlashDeconVsRealPerSpectrum
    {
        private const string Mzml = @"E:\Projects\JurkatTopDown\02-18-20_jurkat_td_rep2_fract10.mzML";
        private const string Msalign = @"E:\Projects\JurkatTopDown\FlashDecon\02-18-20_jurkat_td_rep2_fract10_ms1.msalign";
        private const string ReportPath = @"E:\CodeReview\MetaFlashDecon\deliverables\perspectrum_vs_real.txt";

        private const double MatchPpm = 15.0;
        private const double IsoDa = 1.0033548;
        private const int SampleEvery = 20; // ~160 scans spread evenly across the whole run (every-scan ≈3200 is impractical)

        // Diff OUR candidate set vs OpenMS's (CAND lines in fd_trace_z1_60.txt) to localize the
        // candidate over-count (target: match OpenMS's 5840 before touching scoring).
        [Test]
        public void Candidates_VsOpenMS()
        {
            Assume.That(File.Exists(Mzml), $"missing {Mzml}");
            string trace = @"E:\CodeReview\MetaFlashDecon\difftest\avg_snip\fd_trace_z1_60.txt";
            Assume.That(File.Exists(trace), $"missing {trace}");
            var oms = File.ReadLines(trace).Where(l => l.StartsWith("CAND"))
                .Select(l => { var t = l.Split('\t'); return (m: double.Parse(t[1], CultureInfo.InvariantCulture), z1: int.Parse(t[2]), z2: int.Parse(t[3])); }).ToList();

            var dataFile = MsDataFileReader.GetDataFile(Mzml).LoadAllStaticData();
            var densest = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1).OrderByDescending(s => s.MassSpectrum.Size).First();
            var p = new MetaFlashDeconParameters(minCharge: 1, maxCharge: 60);
            var avg = MetaFlashDeconAveragine.For(p.AverageResidueModel, IsoDa);
            var logPeaks = MetaFlashDeconAlgorithm.BuildLogMzPeaks(densest.MassSpectrum, densest.MassSpectrum.Range, Polarity.Positive);
            double binMul = p.TolDivFactor / (p.DeconvolutionTolerancePpm * 1e-6);
            var ours = MetaFlashDeconCandidateFinder.FindCandidates(logPeaks, 60, new[] { 2, 3, 5, 7, 11 }, binMul, p, avg);

            var omsMono = oms.Select(o => o.m).OrderBy(x => x).ToList();
            var ourMono = ours.Select(c => c.Mass).OrderBy(x => x).ToList();
            int OmsDistinct = oms.Select(o => Math.Round(o.m, 1)).Distinct().Count();
            int OurDistinct = ours.Select(c => Math.Round(c.Mass, 1)).Distinct().Count();
            // how many of ours land within 10ppm of an OpenMS candidate (two-pointer over sorted)
            int matched = 0, jp = 0;
            foreach (double m in ourMono)
            {
                double tolDa = m * 10e-6;
                while (jp < omsMono.Count && omsMono[jp] < m - tolDa) jp++;
                if (jp < omsMono.Count && Math.Abs(omsMono[jp] - m) <= tolDa) matched++;
            }
            TestContext.Progress.WriteLine($"OpenMS cand={oms.Count} (distinct0.1={OmsDistinct})  ours={ours.Count} (distinct0.1={OurDistinct})  ours-within-10ppm-of-OMS={matched}  ours-extra={ours.Count - matched}");
            // mass histogram of ALL ours vs OpenMS (2 kDa buckets)
            string Hist(IEnumerable<double> ms) => string.Join(" ", ms.GroupBy(m => (int)(m / 2000)).OrderBy(g => g.Key).Select(g => $"{g.Key * 2}k={g.Count()}"));
            TestContext.Progress.WriteLine("OpenMS by-2kDa: " + Hist(omsMono));
            TestContext.Progress.WriteLine("OURS   by-2kDa: " + Hist(ourMono));
        }

        // For each of OpenMS's 17 ground-truth masses (fd_cpp_z1_60.txt from flashdeconv_snip),
        // report our nearest CANDIDATE (does candidate gen form it?) and our nearest FINAL group
        // (does scoring keep it?). Pinpoints candidate-gen vs scoring as the bug.
        [Test]
        public void DenseScan_OursVsOpenMSTruth()
        {
            Assume.That(File.Exists(Mzml), $"missing {Mzml}");
            string fdOut = @"E:\CodeReview\MetaFlashDecon\difftest\avg_snip\fd_cpp_z1_60.txt";
            Assume.That(File.Exists(fdOut), $"missing {fdOut}");
            var omsMasses = File.ReadLines(fdOut).Where(l => !l.StartsWith("#") && l.Trim().Length > 0)
                .Select(l => { var t = l.Split('\t'); return (m: double.Parse(t[0], CultureInfo.InvariantCulture),
                    z1: int.Parse(t[1]), z2: int.Parse(t[2])); }).OrderBy(x => x.m).ToList();

            var dataFile = MsDataFileReader.GetDataFile(Mzml).LoadAllStaticData();
            var densest = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1)
                .OrderByDescending(s => s.MassSpectrum.Size).First();
            var p = new MetaFlashDeconParameters(minCharge: 1, maxCharge: 60);
            var avg = MetaFlashDeconAveragine.For(p.AverageResidueModel, IsoDa);
            var logPeaks = MetaFlashDeconAlgorithm.BuildLogMzPeaks(densest.MassSpectrum, densest.MassSpectrum.Range, Polarity.Positive);
            double binMul = p.TolDivFactor / (p.DeconvolutionTolerancePpm * 1e-6);
            var candidates = MetaFlashDeconCandidateFinder.FindCandidates(logPeaks, 60, new[] { 2, 3, 5, 7, 11 }, binMul, p, avg);
            var algo = new MetaFlashDeconAlgorithm(p);
            var scored = algo.ScoreCandidatesViaPeakGroups(densest.MassSpectrum, candidates, p, avg);
            var afterCE = MetaFlashDeconPeakGroup.RemoveChargeErrorPeakGroups(scored, Polarity.Positive);
            double overlapWindow = p.DeconvolutionTolerancePpm * 1e-6 * p.TolDivFactor * p.OverlapDedupTolFactor;
            var final = MetaFlashDeconPeakGroup.RemoveOverlappingPeakGroups(afterCE, overlapWindow);

            TestContext.Progress.WriteLine($"OpenMS truth={omsMasses.Count}  ourCandidates={candidates.Count}  ourFinal={final.Count}");
            foreach (var t in omsMasses)
            {
                var cand = candidates.Where(c => Math.Abs(c.Mass - t.m) <= 2.5).OrderBy(c => Math.Abs(c.Mass - t.m)).FirstOrDefault();
                var fin = final.Where(g => Math.Abs(g.MonoisotopicMass - t.m) <= 2.5).OrderBy(g => Math.Abs(g.MonoisotopicMass - t.m)).FirstOrDefault();
                string candStr = cand != null ? $"{cand.Mass:F1} z{cand.MinAbsCharge}-{cand.MaxAbsCharge}" : "NONE";
                string finStr = fin != null ? $"{fin.MonoisotopicMass:F1} z{fin.MinAbsCharge}-{fin.MaxAbsCharge} cos{fin.IsotopeCosineScore:F3}" : "NONE";
                TestContext.Progress.WriteLine($"  OMS {t.m,9:F1} z{t.z1}-{t.z2} | ourCand {candStr,-22} | ourFinal {finStr}");
            }
            TestContext.Progress.WriteLine("--- ALL our final masses (sorted) ---");
            foreach (var g in final.OrderBy(g => g.MonoisotopicMass))
                TestContext.Progress.WriteLine($"  {g.MonoisotopicMass,9:F1} z{g.MinAbsCharge}-{g.MaxAbsCharge} rep{g.RepAbsCharge} snr={g.Snr:F2} cos={g.IsotopeCosineScore:F3} q={g.QscoreValue:F3}");
            TestContext.Progress.WriteLine($"--- ALL our SCORED (post-gates, pre-CE/overlap, n={scored.Count}) ---");
            foreach (var g in scored.OrderBy(g => g.MonoisotopicMass))
                TestContext.Progress.WriteLine($"  {g.MonoisotopicMass,9:F1} z{g.MinAbsCharge}-{g.MaxAbsCharge} rep{g.RepAbsCharge} snr={g.Snr:F2} cos={g.IsotopeCosineScore:F3} q={g.QscoreValue:F3}");
            TestContext.Progress.WriteLine($"--- ALL our afterCE (n={afterCE.Count}) ---");
            foreach (var g in afterCE.OrderBy(g => g.MonoisotopicMass))
                TestContext.Progress.WriteLine($"  {g.MonoisotopicMass,9:F1} z{g.MinAbsCharge}-{g.MaxAbsCharge} rep{g.RepAbsCharge} snr={g.Snr:F2}");
        }

        // Differential test of OUR MetaFlashDeconAveragine vs OpenMS PrecalculatedAveragine
        // (linked extract: difftest\avg_snip\build\Release\avg_snip.exe -> avg_cpp_out.txt).
        // Compares apex / leftCountFromApex / rightCountFromApex / averageMassDelta / the
        // normalised b vector, per mass, to confirm/refute the averagine boundary.
        // Writes the densest MS1 scan's centroid peaks ("mz intensity" per line) for the
        // flashdeconv_snip OpenMS extract to consume.
        [Test]
        public void DumpDensestScanPeaks()
        {
            Assume.That(File.Exists(Mzml), $"missing {Mzml}");
            var dataFile = MsDataFileReader.GetDataFile(Mzml).LoadAllStaticData();
            var densest = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1)
                .OrderByDescending(s => s.MassSpectrum.Size).First();
            var sb = new System.Text.StringBuilder();
            var x = densest.MassSpectrum.XArray; var y = densest.MassSpectrum.YArray;
            for (int i = 0; i < x.Length; i++)
                sb.Append(x[i].ToString("R", CultureInfo.InvariantCulture)).Append(' ')
                  .Append(y[i].ToString("R", CultureInfo.InvariantCulture)).Append('\n');
            string outp = @"E:\CodeReview\MetaFlashDecon\difftest\avg_snip\scan2301_peaks.txt";
            File.WriteAllText(outp, sb.ToString());
            TestContext.Progress.WriteLine($"wrote {densest.MassSpectrum.Size} peaks of scan {densest.OneBasedScanNumber} -> {outp}");
        }

        // Scans to fidelity-check vs the REAL OpenMS (flashdeconv_snip, our exact z1-60 params):
        // sparse / medium / dense, picked from the per-spectrum report.
        private static readonly int[] SnipScans = { 1332, 2625, 3551 };
        private static string SnipDir => @"E:\CodeReview\MetaFlashDecon\difftest\avg_snip";

        [Test]
        public void DumpScansForSnip()
        {
            Assume.That(File.Exists(Mzml), $"missing {Mzml}");
            var dataFile = MsDataFileReader.GetDataFile(Mzml).LoadAllStaticData();
            foreach (int scanNum in SnipScans)
            {
                var scan = dataFile.GetAllScansList().First(s => s.OneBasedScanNumber == scanNum);
                var sb = new System.Text.StringBuilder();
                var x = scan.MassSpectrum.XArray; var y = scan.MassSpectrum.YArray;
                for (int i = 0; i < x.Length; i++)
                    sb.Append(x[i].ToString("R", CultureInfo.InvariantCulture)).Append(' ')
                      .Append(y[i].ToString("R", CultureInfo.InvariantCulture)).Append('\n');
                File.WriteAllText($@"{SnipDir}\scan{scanNum}_peaks.txt", sb.ToString());
                TestContext.Progress.WriteLine($"wrote {scan.MassSpectrum.Size} peaks of scan {scanNum}");
            }
            TestContext.Progress.WriteLine("Now run: flashdeconv_snip scan{N}_peaks.txt fd_{N}.txt 1 60  for each, then CompareScansVsSnip.");
        }

        [Test]
        public void CompareScansVsSnip()
        {
            Assume.That(File.Exists(Mzml), $"missing {Mzml}");
            var p = new MetaFlashDeconParameters(minCharge: 1, maxCharge: 60);
            var dataFile = MsDataFileReader.GetDataFile(Mzml).LoadAllStaticData();
            foreach (int scanNum in SnipScans)
            {
                string fd = $@"{SnipDir}\fd_{scanNum}.txt";
                if (!File.Exists(fd)) { TestContext.Progress.WriteLine($"scan {scanNum}: missing {fd} (run flashdeconv_snip)"); continue; }
                var oms = File.ReadLines(fd).Where(l => l.Trim().Length > 0 && !l.StartsWith("#"))
                    .Select(l => double.Parse(l.Split('\t')[0], CultureInfo.InvariantCulture)).OrderBy(m => m).ToList();
                var scan = dataFile.GetAllScansList().First(s => s.OneBasedScanNumber == scanNum);
                var ours = Deconvoluter.Deconvolute(scan.MassSpectrum, p).Select(e => e.MonoisotopicMass).OrderBy(m => m).ToList();

                var omsUsed = new bool[oms.Count];
                int matched = 0;
                foreach (double m in ours)
                {
                    double tol = m * 10e-6; int best = -1; double bestErr = double.MaxValue;
                    for (int i = 0; i < oms.Count; i++)
                    { if (omsUsed[i]) continue; double d = Math.Abs(oms[i] - m); if (d <= tol && d < bestErr) { bestErr = d; best = i; } }
                    if (best >= 0) { omsUsed[best] = true; matched++; }
                }
                int extra = ours.Count - matched, missed = omsUsed.Count(u => !u);
                TestContext.Progress.WriteLine($"scan {scanNum}: openms={oms.Count} ours={ours.Count} | matched={matched} extra={extra} missed={missed}");
            }
        }

        [Test]
        public void DiagnoseScan1332Misses()
        {
            int scanNum = 1332;
            string fd = $@"{SnipDir}\fd_{scanNum}.txt";
            Assume.That(File.Exists(Mzml) && File.Exists(fd), $"need {Mzml} + {fd}");
            var oms = File.ReadLines(fd).Where(l => l.Trim().Length > 0 && !l.StartsWith("#"))
                .Select(l => { var t = l.Split('\t'); return (m: double.Parse(t[0], CultureInfo.InvariantCulture),
                    z1: t.Length > 1 ? int.Parse(t[1]) : 0, z2: t.Length > 2 ? int.Parse(t[2]) : 0); }).OrderBy(x => x.m).ToList();

            var dataFile = MsDataFileReader.GetDataFile(Mzml).LoadAllStaticData();
            var scan = dataFile.GetAllScansList().First(s => s.OneBasedScanNumber == scanNum);
            var p = new MetaFlashDeconParameters(minCharge: 1, maxCharge: 60);
            var avg = MetaFlashDeconAveragine.For(p.AverageResidueModel, MetaFlashDeconAlgorithm.IsoDaDistance55K);
            var logPeaks = MetaFlashDeconAlgorithm.BuildLogMzPeaks(scan.MassSpectrum, scan.MassSpectrum.Range, Polarity.Positive);
            double binMul = p.TolDivFactor / (p.DeconvolutionTolerancePpm * 1e-6);
            var candidates = MetaFlashDeconCandidateFinder.FindCandidates(logPeaks, 60, new[] { 2, 3, 5, 7, 11 }, binMul, p, avg);
            var algo = new MetaFlashDeconAlgorithm(p);
            var scored = algo.ScoreCandidatesViaPeakGroups(scan.MassSpectrum, candidates, p, avg);
            var afterCE = MetaFlashDeconPeakGroup.RemoveChargeErrorPeakGroups(scored, Polarity.Positive);
            double ow = p.DeconvolutionTolerancePpm * 1e-6 * p.TolDivFactor * p.OverlapDedupTolFactor;
            var final = MetaFlashDeconPeakGroup.RemoveOverlappingPeakGroups(afterCE, ow);

            TestContext.Progress.WriteLine($"scan {scanNum}: openms={oms.Count} cands={candidates.Count} scored={scored.Count} afterCE={afterCE.Count} final={final.Count}");
            string N(IEnumerable<double> ms, double m) { var h = ms.Where(x => Math.Abs(x - m) <= 2.5).OrderBy(x => Math.Abs(x - m)).Select(x => (double?)x).FirstOrDefault(); return h.HasValue ? h.Value.ToString("F1") : "—"; }
            foreach (var t in oms)
            {
                string c = N(candidates.Select(x => x.Mass), t.m);
                string s = N(scored.Select(x => x.MonoisotopicMass), t.m);
                string ce = N(afterCE.Select(x => x.MonoisotopicMass), t.m);
                string f = N(final.Select(x => x.MonoisotopicMass), t.m);
                bool found = final.Any(x => Math.Abs(x.MonoisotopicMass - t.m) <= t.m * 10e-6);
                TestContext.Progress.WriteLine($"  {(found ? " " : "X")} OMS {t.m,9:F1} z{t.z1}-{t.z2} | cand {c,-9} scored {s,-9} afterCE {ce,-9} final {f}");
            }
        }

        [Test]
        public void Averagine_VsOpenMS()
        {
            string cppOut = @"E:\CodeReview\MetaFlashDecon\difftest\avg_snip\avg_cpp_out.txt";
            Assume.That(File.Exists(cppOut), $"missing {cppOut} (run avg_snip.exe)");

            var avg = MetaFlashDeconAveragine.For(new Averagine(), IsoDa);
            foreach (var raw in File.ReadLines(cppOut))
            {
                var t = raw.Split('\t');
                double m = double.Parse(t[1], CultureInfo.InvariantCulture);
                int cApex = int.Parse(t[3]), cLeft = int.Parse(t[5]), cRight = int.Parse(t[7]);
                double cAvgDelta = double.Parse(t[9], CultureInfo.InvariantCulture);
                int bIdx = Array.IndexOf(t, "b");
                var cB = t.Skip(bIdx + 1).Select(x => double.Parse(x, CultureInfo.InvariantCulture)).ToArray();

                int oApex = avg.GetApexIndex(m), oLeft = avg.GetLeftCountFromApex(m), oRight = avg.GetRightCountFromApex(m);
                double oAvgDelta = avg.GetAverageMassDelta(m);
                var oB = avg.Get(m);

                // cosine + max abs diff of the b vectors (aligned at index 0 = monoisotopic)
                int n = Math.Max(cB.Length, oB.Length);
                double dot = 0, na = 0, nb = 0, maxd = 0;
                for (int k = 0; k < n; k++)
                {
                    double a = k < oB.Length ? oB[k] : 0, b = k < cB.Length ? cB[k] : 0;
                    dot += a * b; na += a * a; nb += b * b; maxd = Math.Max(maxd, Math.Abs(a - b));
                }
                double cos = (na > 0 && nb > 0) ? dot / Math.Sqrt(na * nb) : 0;
                TestContext.Progress.WriteLine(
                    $"m={m,8:F1}  apex {oApex}/{cApex}  left {oLeft}/{cLeft}  right {oRight}/{cRight}  " +
                    $"avgDelta {oAvgDelta:F3}/{cAvgDelta:F3}  n {oB.Length}/{cB.Length}  Bcos={cos:F4} maxBdiff={maxd:F4}  (ours/openms)");
            }
        }

        [Test]
        public void DenseScan_ScoringStageCounts()
        {
            Assume.That(File.Exists(Mzml), $"missing {Mzml}");
            var dataFile = MsDataFileReader.GetDataFile(Mzml).LoadAllStaticData();
            var densest = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1)
                .OrderByDescending(s => s.MassSpectrum.Size).First();

            var p = new MetaFlashDeconParameters(minCharge: 1, maxCharge: 60);
            var avg = MetaFlashDeconAveragine.For(p.AverageResidueModel, IsoDa);
            var logPeaks = MetaFlashDeconAlgorithm.BuildLogMzPeaks(densest.MassSpectrum, densest.MassSpectrum.Range, Polarity.Positive);
            double binMul = p.TolDivFactor / (p.DeconvolutionTolerancePpm * 1e-6);
            int[] harmonicCharges = { 2, 3, 5, 7, 11 };

            var candidates = MetaFlashDeconCandidateFinder.FindCandidates(logPeaks, 60, harmonicCharges, binMul, p, avg);
            var algo = new MetaFlashDeconAlgorithm(p);
            var scored = algo.ScoreCandidatesViaPeakGroups(densest.MassSpectrum, candidates, p, avg);
            var afterCE = MetaFlashDeconPeakGroup.RemoveChargeErrorPeakGroups(scored, Polarity.Positive);
            double overlapWindow = p.DeconvolutionTolerancePpm * 1e-6 * p.TolDivFactor * p.OverlapDedupTolFactor;
            var final = MetaFlashDeconPeakGroup.RemoveOverlappingPeakGroups(afterCE, overlapWindow);

            // distinct masses among the final (to detect un-merged near-duplicates)
            var distinct = final.Select(g => Math.Round(g.MonoisotopicMass, 1)).Distinct().Count();
            TestContext.Progress.WriteLine(
                $"densest scan {densest.OneBasedScanNumber} ({densest.MassSpectrum.Size} peaks): " +
                $"candidates={candidates.Count} -> afterScoringGates={scored.Count} -> afterChargeError={afterCE.Count} -> afterOverlap={final.Count} (distinct masses={distinct})");
        }

        [Test]
        public void DenseScan_ClassifyExtrasVsReal()
        {
            Assume.That(File.Exists(Mzml), $"missing {Mzml}");
            Assume.That(File.Exists(Msalign), $"missing {Msalign}");
            var dataFile = MsDataFileReader.GetDataFile(Mzml).LoadAllStaticData();
            var densest = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1)
                .OrderByDescending(s => s.MassSpectrum.Size).First();
            int scanNo = densest.OneBasedScanNumber;
            var real = ParseMsalign(Msalign, out _);
            var realMasses = real.TryGetValue(scanNo, out var rm) ? rm.OrderBy(m => m).ToList() : new List<double>();

            var p = new MetaFlashDeconParameters(minCharge: 1, maxCharge: 60);
            var avg = MetaFlashDeconAveragine.For(p.AverageResidueModel, IsoDa);
            var logPeaks = MetaFlashDeconAlgorithm.BuildLogMzPeaks(densest.MassSpectrum, densest.MassSpectrum.Range, Polarity.Positive);
            double binMul = p.TolDivFactor / (p.DeconvolutionTolerancePpm * 1e-6);
            var candidates = MetaFlashDeconCandidateFinder.FindCandidates(logPeaks, 60, new[] { 2, 3, 5, 7, 11 }, binMul, p, avg);
            var candWidth = candidates.GroupBy(c => Math.Min(c.MaxAbsCharge - c.MinAbsCharge, 5))
                .OrderBy(g => g.Key).Select(g => $"span{g.Key}={g.Count()}");
            TestContext.Progress.WriteLine($"candidate charge-span histogram (n={candidates.Count}): " + string.Join(" ", candWidth));
            var algo = new MetaFlashDeconAlgorithm(p);
            var scored = algo.ScoreCandidatesViaPeakGroups(densest.MassSpectrum, candidates, p, avg);
            // ── instrument removeChargeError peak-sharing on the scored (pre-CE) groups ──
            {
                double cm = Chemistry.Constants.ProtonMass;
                var p2g = new Dictionary<double, HashSet<int>>();
                var mzi = new Dictionary<double, double>();
                for (int i = 0; i < scored.Count; i++)
                    foreach (var pk in scored[i].SignalPeaks)
                    { if (!p2g.TryGetValue(pk.Mz, out var s)) { s = new HashSet<int>(); p2g[pk.Mz] = s; } s.Add(i); mzi[pk.Mz] = pk.Intensity; }
                int sharedPeaks = p2g.Count(kv => kv.Value.Count > 1);
                var ovl = new double[scored.Count];
                foreach (var kv in p2g)
                {
                    if (kv.Value.Count == 1) continue;
                    double pmz = kv.Key, pint = mzi[pmz];
                    foreach (int i in kv.Value)
                    {
                        bool isOv = false; int rz1 = (int)Math.Round(scored[i].MonoisotopicMass / (pmz - cm), MidpointRounding.AwayFromZero);
                        foreach (int j in kv.Value)
                        {
                            if (i == j) continue;
                            int rz2 = (int)Math.Round(scored[j].MonoisotopicMass / (pmz - cm), MidpointRounding.AwayFromZero);
                            if (rz1 == rz2) continue;
                            if (scored[i].GetChargeSnr(rz1) > scored[j].GetChargeSnr(rz2) * 2.0) continue;
                            isOv = true; break;
                        }
                        if (isOv) ovl[i] += pint;
                    }
                }
                int groupsSharing = Enumerable.Range(0, scored.Count).Count(i => scored[i].SignalPeaks.Any(pk => p2g[pk.Mz].Count > 1));
                int wouldDrop = Enumerable.Range(0, scored.Count).Count(i => ovl[i] >= scored[i].GetIntensity() * 0.5);
                var ratios = Enumerable.Range(0, scored.Count).Select(i => scored[i].GetIntensity() > 0 ? ovl[i] / scored[i].GetIntensity() : 0).OrderBy(r => r).ToList();
                TestContext.Progress.WriteLine($"CE-instrument: scored={scored.Count} distinctPeaks={p2g.Count} sharedPeaks={sharedPeaks} groupsSharingAnyPeak={groupsSharing} wouldDrop(>=0.5)={wouldDrop}");
                TestContext.Progress.WriteLine($"  overlap-ratio pct: p10={ratios[ratios.Count/10]:F2} p50={ratios[ratios.Count/2]:F2} p90={ratios[ratios.Count*9/10]:F2} max={ratios.Last():F2}");
                // Dump 27889.7's z32 signal peaks (to compare against OpenMS - does it include 872.84?)
                {
                    int t27 = -1;
                    for (int i = 0; i < scored.Count; i++)
                        if (Math.Abs(scored[i].MonoisotopicMass - 27889.7) < 0.5) { t27 = i; break; }
                    if (t27 >= 0)
                    {
                        double monoExact = scored[t27].MonoisotopicMass;
                        double cmzExact = monoExact / 32 + Chemistry.Constants.ProtonMass;
                        double isoDeltaExact = MetaFlashDeconAlgorithm.IsoDaDistance55K / 32;
                        TestContext.Progress.WriteLine($"OUR 27889.7 exact: mono={monoExact:F8} cmz_z32={cmzExact:F8} iso_delta_z32={isoDeltaExact:F10}");
                        // For peak 873.348 (iso 25): compute exact signal check
                        double pmzCheck = 873.3480;
                        int isoCheck = (int)Math.Round((pmzCheck - cmzExact) / isoDeltaExact, MidpointRounding.AwayFromZero);
                        double errCheck = Math.Abs(pmzCheck - cmzExact - isoCheck * isoDeltaExact);
                        double tolCheck = pmzCheck * 4e-6;
                        TestContext.Progress.WriteLine($"  Signal-check at pmz=873.348: iso={isoCheck} error={errCheck:F6} threshold={tolCheck:F6} signal={errCheck <= tolCheck}");
                        TestContext.Progress.WriteLine($"OUR 27889.7 z32 signal peaks:");
                        foreach (var pk in scored[t27].SignalPeaks.Where(p => p.AbsCharge == 32).OrderBy(p => p.Mz))
                            TestContext.Progress.WriteLine($"  z32 mz={pk.Mz:F4} int={pk.Intensity:E2} iso={pk.IsotopeIndex}");
                    }
                }
                // Per-charge SNRs for 27889.7 (compare against OpenMS SCORED dump)
                {
                    int t27 = -1;
                    for (int i = 0; i < scored.Count; i++)
                        if (Math.Abs(scored[i].MonoisotopicMass - 27889.7) < 0.5) { t27 = i; break; }
                    if (t27 >= 0)
                    {
                        var sb = new System.Text.StringBuilder($"OUR 27889.7 per-charge SNRs: ");
                        for (int z = scored[t27].MinAbsCharge; z <= scored[t27].MaxAbsCharge; z++)
                            sb.Append($"z{z}:{scored[t27].GetChargeSnr(z):F4} ");
                        TestContext.Progress.WriteLine(sb.ToString());
                        sb.Clear(); sb.Append($"OUR 6970.7 per-charge SNRs: ");
                        int t697 = -1;
                        for (int i = 0; i < scored.Count; i++)
                            if (Math.Abs(scored[i].MonoisotopicMass - 6970.7) < 0.5) { t697 = i; break; }
                        if (t697 >= 0)
                            for (int z = scored[t697].MinAbsCharge; z <= scored[t697].MaxAbsCharge; z++)
                                sb.Append($"z{z}:{scored[t697].GetChargeSnr(z):F4} ");
                        TestContext.Progress.WriteLine(sb.ToString());
                    }
                }
                // Trace contributions to 6970.7's overlap (the first dup index)
                int tgt = -1;
                for (int i = 0; i < scored.Count; i++)
                    if (Math.Abs(scored[i].MonoisotopicMass - 6970.7) < 0.5) { tgt = i; break; }
                if (tgt >= 0)
                {
                    TestContext.Progress.WriteLine($"  TRACE contributions to 6970.7 (idx={tgt}, mass={scored[tgt].MonoisotopicMass:F4}):");
                    foreach (var kv in p2g)
                    {
                        if (!kv.Value.Contains(tgt) || kv.Value.Count == 1) continue;
                        double pmz2 = kv.Key, pint2 = mzi[pmz2];
                        int rz1 = (int)Math.Round(scored[tgt].MonoisotopicMass / (pmz2 - cm), MidpointRounding.AwayFromZero);
                        foreach (int j in kv.Value)
                        {
                            if (j == tgt) continue;
                            int rz2 = (int)Math.Round(scored[j].MonoisotopicMass / (pmz2 - cm), MidpointRounding.AwayFromZero);
                            string match = rz1 == rz2 ? "rzMATCH" : "rzDIFF";
                            string snrCmp = rz1 == rz2 ? "" : ($"snr_i({rz1})={scored[tgt].GetChargeSnr(rz1):F2} snr_j({rz2})={scored[j].GetChargeSnr(rz2):F2} i>2j? {(scored[tgt].GetChargeSnr(rz1) > 2 * scored[j].GetChargeSnr(rz2))}");
                            if (match == "rzDIFF") TestContext.Progress.WriteLine($"    pmz={pmz2:F4} pint={pint2:E2} vs j={j}(mass={scored[j].MonoisotopicMass:F2}) {match} {snrCmp}");
                        }
                    }
                }
                TestContext.Progress.WriteLine("  PER-GROUP (mono z rep ovl intensity ratio dropped):");
                for (int i = 0; i < scored.Count; i++)
                {
                    double inten = scored[i].GetIntensity();
                    double ratio = inten > 0 ? ovl[i] / inten : 0;
                    string mark = ratio >= 0.5 ? "DROP" : "keep";
                    TestContext.Progress.WriteLine($"    {scored[i].MonoisotopicMass,9:F1} z{scored[i].MinAbsCharge}-{scored[i].MaxAbsCharge} rep{scored[i].RepAbsCharge} ovl={ovl[i]:E2} int={inten:E2} ratio={ratio:F3} {mark}");
                }
            }
            var afterCE = MetaFlashDeconPeakGroup.RemoveChargeErrorPeakGroups(scored, Polarity.Positive);
            double overlapWindow = p.DeconvolutionTolerancePpm * 1e-6 * p.TolDivFactor * p.OverlapDedupTolFactor;
            var final = MetaFlashDeconPeakGroup.RemoveOverlappingPeakGroups(afterCE, overlapWindow).OrderBy(g => g.MonoisotopicMass).ToList();

            // classify each surviving group vs the real masses for this scan
            string Classify(double m)
            {
                foreach (var r in realMasses) if (Math.Abs(m - r) <= r * 15e-6) return "matched";
                foreach (var r in realMasses) for (int k = 1; k <= 3; k++)
                    if (Math.Abs(m - (r + k * IsoDa)) <= r * 15e-6 || Math.Abs(m - (r - k * IsoDa)) <= r * 15e-6) return "iso-off";
                foreach (var r in realMasses) for (int k = 2; k <= 6; k++)
                    if (Math.Abs(m - r * k) <= m * 15e-6 || Math.Abs(m - r / k) <= m * 15e-6) return "harmonic";
                return "spurious";
            }
            var groups = final.GroupBy(g => Classify(g.MonoisotopicMass)).ToDictionary(g => g.Key, g => g.ToList());

            void Stat(string name, Func<MetaFlashDeconPeakGroup, double> sel)
            {
                foreach (var cls in new[] { "matched", "iso-off", "harmonic", "spurious" })
                {
                    if (!groups.TryGetValue(cls, out var gs) || gs.Count == 0) continue;
                    var vals = gs.Select(sel).OrderBy(v => v).ToList();
                    TestContext.Progress.WriteLine($"  {name,-10} {cls,-9} n={vals.Count,3}  min={vals.First():F3} med={vals[vals.Count/2]:F3} max={vals.Last():F3}");
                }
            }
            TestContext.Progress.WriteLine($"scan {scanNo}: real masses={realMasses.Count}, our survivors={final.Count}");
            TestContext.Progress.WriteLine("  class counts: " + string.Join("  ", groups.OrderByDescending(g=>g.Value.Count).Select(g => $"{g.Key}={g.Value.Count}")));
            Stat("cosine", g => g.IsotopeCosineScore);
            Stat("snr", g => g.Snr);
            Stat("qscore", g => g.QscoreValue);
            Stat("chargefit", g => g.ChargeScore);
            Stat("chargeN", g => g.MaxAbsCharge - g.MinAbsCharge + 1);
            Stat("mass", g => g.MonoisotopicMass);

            TestContext.Progress.WriteLine("REAL masses (" + realMasses.Count + "): " +
                string.Join(" ", realMasses.Select(m => m.ToString("F1"))));
            TestContext.Progress.WriteLine("OUR masses (mono@repZ, snr) sorted:");
            foreach (var g in final)
                TestContext.Progress.WriteLine($"  {g.MonoisotopicMass,10:F1} z{g.MinAbsCharge}-{g.MaxAbsCharge} rep{g.RepAbsCharge} snr={g.Snr:F2} cos={g.IsotopeCosineScore:F3} int={g.GetIntensity():E2}");
        }

        [Test]
        public void ComparePerSpectrumToRealFlashDeconv()
        {
            Assume.That(File.Exists(Mzml), $"missing {Mzml}");
            Assume.That(File.Exists(Msalign), $"missing {Msalign}");

            var real = ParseMsalign(Msalign, out var realRtSec); // scan -> masses; + scan -> RT(seconds)

            var dataFile = MsDataFileReader.GetDataFile(Mzml).LoadAllStaticData();
            var ms1 = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1).ToDictionary(s => s.OneBasedScanNumber);

            var p = new MetaFlashDeconParameters(minCharge: 1, maxCharge: 60);

            var log = new List<string>();
            void Log(string s) { TestContext.Progress.WriteLine(s); log.Add(s); }

            Log($"per-spectrum MetaFlashDecon vs real FLASHDeconv  ({DateTime.Now:HH:mm:ss})");
            Log($"match within {MatchPpm} ppm; ±1 Da = isotope-offset error. params minZ=1 maxZ=60.");
            Log("scan | realN ourN | matched +-1Da extra missed | notes");

            var scans = real.Keys.Where(k => ms1.ContainsKey(k)).OrderBy(k => k).ToList();
            int tReal = 0, tOur = 0, tMatch = 0, tPlus1 = 0, tExtra = 0, tMissed = 0, nScans = 0;
            int tLowReal = 0, tLowOur = 0; // < 5000 Da

            for (int si = 0; si < scans.Count; si += SampleEvery)
            {
                int scan = scans[si];
                var realMasses = real[scan];
                var ours = Deconvoluter.Deconvolute(ms1[scan].MassSpectrum, p)
                    .Select(e => e.MonoisotopicMass).OrderBy(m => m).ToList();

                int matched = 0, plus1 = 0, extra = 0;
                var realUsed = new bool[realMasses.Count];
                foreach (double om in ours)
                {
                    int hit = FindMatch(realMasses, om, MatchPpm, realUsed, exactOnly: true);
                    if (hit >= 0) { matched++; realUsed[hit] = true; continue; }
                    int hit1 = FindMatch(realMasses, om, MatchPpm, realUsed, exactOnly: false);
                    if (hit1 >= 0) { plus1++; realUsed[hit1] = true; continue; }
                    extra++;
                }
                int missed = realUsed.Count(u => !u);

                tReal += realMasses.Count; tOur += ours.Count; tMatch += matched;
                tPlus1 += plus1; tExtra += extra; tMissed += missed; nScans++;
                tLowReal += realMasses.Count(m => m < 5000);
                tLowOur += ours.Count(m => m < 5000);

                double ourRtMin = ms1[scan].RetentionTime;
                double realRtMin = realRtSec.TryGetValue(scan, out var rs) ? rs / 60.0 : -1;
                if (si % (SampleEvery * 5) == 0 || ours.Count > realMasses.Count * 2)
                    Log($"{scan,5} | {realMasses.Count,5} {ours.Count,4} | {matched,7} {plus1,5} {extra,5} {missed,6} | ourRT={ourRtMin:F2} realRT={realRtMin:F2}");
            }

            Log("");
            Log($"TOT: scans={nScans}  realMasses={tReal}  ourMasses={tOur}  (our/real={(double)tOur/Math.Max(1,tReal):F2}x)");
            Log($"     matched={tMatch}  +-1Da-off={tPlus1}  our-extra={tExtra}  real-missed={tMissed}");
            Log($"     <5kDa: real={tLowReal} ours={tLowOur}   (low-mass recovery)");
            File.WriteAllText(ReportPath, string.Join(Environment.NewLine, log) + Environment.NewLine);
            Log($"wrote {ReportPath}");
        }

        // ── Feature-tracer "piece A" differential test (vs real OpenMS MassFeatureTrace math) ──
        // Runs our decon on a contiguous block of real MS1 scans, clusters the surviving peak groups
        // into neutral-mass traces, and for each trace computes the per-isotope-accumulation +
        // averagine-cosine + isotope-offset EXACTLY as OpenMS MassFeatureTrace::findFeatures does
        // (cpp:84-149). Writes the trace inputs (mft_traces.txt) for massfeature_snip and our results
        // (mft_ours.txt). Run massfeature_snip on mft_traces.txt, then FeatureCosine_VsOpenMS compares.
        private const string MftTraces = @"E:\CodeReview\MetaFlashDecon\difftest\avg_snip\mft_traces.txt";
        private const string MftOurs   = @"E:\CodeReview\MetaFlashDecon\difftest\avg_snip\mft_ours.txt";
        private const string MftCpp    = @"E:\CodeReview\MetaFlashDecon\difftest\avg_snip\mft_cpp_out.txt";

        [Test]
        public void FeatureCosine_DumpForOpenMS()
        {
            Assume.That(File.Exists(Mzml), $"missing {Mzml}");
            double isoDa = MetaFlashDeconAlgorithm.IsoDaDistance55K; // 1.002371 — same as OpenMS pg.getIsotopeDaDistance()
            var p = new MetaFlashDeconParameters(minCharge: 1, maxCharge: 60);
            var avg = MetaFlashDeconAveragine.For(p.AverageResidueModel, isoDa);
            double binMul = p.TolDivFactor / (p.DeconvolutionTolerancePpm * 1e-6);
            int[] harm = { 2, 3, 5, 7, 11 };
            var algo = new MetaFlashDeconAlgorithm(p);

            var dataFile = MsDataFileReader.GetDataFile(Mzml).LoadAllStaticData();
            var ms1 = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1).OrderBy(s => s.OneBasedScanNumber).ToList();
            // contiguous busy block (mid-run), ~150 scans -> real co-eluting traces
            var slice = ms1.Skip(Math.Max(0, ms1.Count / 2 - 75)).Take(150).ToList();

            // per-scan surviving peak groups (carry PerIsotopeInt)
            var all = new List<MetaFlashDeconPeakGroup>();
            foreach (var s in slice)
            {
                var logPeaks = MetaFlashDeconAlgorithm.BuildLogMzPeaks(s.MassSpectrum, s.MassSpectrum.Range, Polarity.Positive);
                var cands = MetaFlashDeconCandidateFinder.FindCandidates(logPeaks, 60, harm, binMul, p, avg);
                var scored = algo.ScoreCandidatesViaPeakGroups(s.MassSpectrum, cands, p, avg);
                var afterCE = MetaFlashDeconPeakGroup.RemoveChargeErrorPeakGroups(scored, Polarity.Positive);
                double ow = p.DeconvolutionTolerancePpm * 1e-6 * p.TolDivFactor * p.OverlapDedupTolFactor;
                all.AddRange(MetaFlashDeconPeakGroup.RemoveOverlappingPeakGroups(afterCE, ow));
            }

            // cluster across scans by mono mass within 15 ppm -> traces (member sets for the A math)
            var traces = new List<List<MetaFlashDeconPeakGroup>>();
            List<MetaFlashDeconPeakGroup> cur = null; double curMass = 0;
            foreach (var g in all.OrderBy(g => g.MonoisotopicMass))
            {
                if (cur != null && Math.Abs(g.MonoisotopicMass - curMass) <= curMass * 15e-6) cur.Add(g);
                else { cur = new List<MetaFlashDeconPeakGroup>(); traces.Add(cur); cur.Add(g); }
                curMass = cur.Average(x => x.MonoisotopicMass);
            }
            var keep = traces.Where(t => t.Count >= 3).Take(400).ToList();

            var sbIn = new System.Text.StringBuilder();
            var sbOurs = new System.Text.StringBuilder();
            int id = 0;
            foreach (var t in keep)
            {
                double totInt = t.Sum(g => g.GetIntensity());
                double centroid = totInt > 0 ? t.Sum(g => g.MonoisotopicMass * g.GetIntensity()) / totInt : t[0].MonoisotopicMass;

                sbIn.Append($"TRACE\t{id}\t{centroid.ToString("R", CultureInfo.InvariantCulture)}\t{isoDa.ToString("R", CultureInfo.InvariantCulture)}\t{t.Count}\n");
                foreach (var g in t)
                {
                    var iso = g.PerIsotopeInt; // double[], min-negative(-1)-based — same as OpenMS getIsotopeIntensities
                    sbIn.Append($"PG\t{g.MonoisotopicMass.ToString("R", CultureInfo.InvariantCulture)}\t{iso.Length}");
                    foreach (var v in iso) sbIn.Append('\t').Append(((float)v).ToString("R", CultureInfo.InvariantCulture)); // float — OpenMS per_isotope is float
                    sbIn.Append('\n');
                }

                // OUR side: replicate MassFeatureTrace.cpp:84-138 accumulation in FLOAT (matches std::vector<float>)
                const int maxIso = 400; // >= snip's avg.getMaxIsotopeIndex(); cosine ignores trailing zeros
                var perIso = new float[maxIso];
                foreach (var g in t)
                {
                    var iso = g.PerIsotopeInt;
                    int isoOff = (int)(0.5 + (g.MonoisotopicMass - centroid) / isoDa); // C# (int) truncates toward 0 == C++ int()
                    for (int i = 0; i < perIso.Length - isoOff; i++)
                    {
                        if (i + isoOff < 0 || i >= iso.Length) continue;
                        perIso[i + isoOff] += (float)iso[i];
                    }
                }
                var perIsoD = new double[maxIso];
                for (int i = 0; i < maxIso; i++) perIsoD[i] = perIso[i];
                double cos = MetaFlashDeconPeakGroup.GetIsotopeCosineAndDetermineIsotopeIndex(
                    perIsoD, avg.Get(centroid), avg.GetApexIndex(centroid), 0, 0, 2, out int off);
                double corrected = centroid + off * isoDa;
                sbOurs.Append($"FEAT\t{id}\t{cos.ToString("R", CultureInfo.InvariantCulture)}\t{off}\t{corrected.ToString("R", CultureInfo.InvariantCulture)}\n");
                id++;
            }
            File.WriteAllText(MftTraces, sbIn.ToString());
            File.WriteAllText(MftOurs, sbOurs.ToString());
            TestContext.Progress.WriteLine($"wrote {keep.Count} traces -> {MftTraces} + {MftOurs}. Now run massfeature_snip, then FeatureCosine_VsOpenMS.");
        }

        [Test]
        public void FeatureCosine_VsOpenMS()
        {
            Assume.That(File.Exists(MftOurs), $"missing {MftOurs} (run FeatureCosine_DumpForOpenMS)");
            Assume.That(File.Exists(MftCpp), $"missing {MftCpp} (run massfeature_snip)");
            (int id, double cos, int off, double mass) Parse(string l)
            { var t = l.Split('\t'); return (int.Parse(t[1]), double.Parse(t[2], CultureInfo.InvariantCulture), int.Parse(t[3]), double.Parse(t[4], CultureInfo.InvariantCulture)); }
            var ours = File.ReadLines(MftOurs).Where(l => l.StartsWith("FEAT")).Select(Parse).ToDictionary(x => x.id);
            var cpp = File.ReadLines(MftCpp).Where(l => l.StartsWith("FEAT")).Select(Parse).ToDictionary(x => x.id);

            double maxCosDiff = 0, maxMassDiff = 0; int offMismatch = 0, crossFilter = 0, n = 0;
            foreach (var kv in cpp)
            {
                if (!ours.TryGetValue(kv.Key, out var o)) continue;
                n++;
                double cosDiff = Math.Abs(o.cos - kv.Value.cos);
                maxCosDiff = Math.Max(maxCosDiff, cosDiff);
                maxMassDiff = Math.Max(maxMassDiff, Math.Abs(o.mass - kv.Value.mass));
                if (o.off != kv.Value.off) offMismatch++;
                if ((o.cos >= 0.75) != (kv.Value.cos >= 0.75)) crossFilter++; // would the 0.75 gate flip?
            }
            TestContext.Progress.WriteLine($"compared n={n}  maxCosDiff={maxCosDiff:E3}  maxMassDiff={maxMassDiff:E3}  offsetMismatch={offMismatch}  cross-0.75-filter={crossFilter}");
            Assert.That(offMismatch, Is.EqualTo(0), "isotope offset diverges -> accumulation/indexing bug");
            Assert.That(maxCosDiff, Is.LessThan(0.02), "isotope-cosine diverges beyond averagine-boundary tolerance");
        }

        // ── MassTraceDetection ("piece C") differential test vs real OpenMS ──
        // Runs our decon on a contiguous block, then DetectTracesDiagnostic returns the per-scan collapsed
        // input peaks + our raw traces. Serialises the input for masstrace_snip (which runs the REAL
        // OpenMS MassTraceDetection on identical peaks) and our traces, then MassTrace_VsOpenMS compares
        // trace-for-trace — isolating port correctness from per-spectrum input differences.
        private const string MtdPeaks = @"E:\CodeReview\MetaFlashDecon\difftest\avg_snip\mtd_peaks.txt";
        private const string MtdOurs  = @"E:\CodeReview\MetaFlashDecon\difftest\avg_snip\mtd_ours.txt";
        private const string MtdCpp   = @"E:\CodeReview\MetaFlashDecon\difftest\avg_snip\mtd_cpp_out.txt";

        [Test]
        public void MassTrace_DumpForOpenMS()
        {
            Assume.That(File.Exists(Mzml), $"missing {Mzml}");
            var p = new MetaFlashDeconParameters(minCharge: 1, maxCharge: 60);
            var dataFile = MsDataFileReader.GetDataFile(Mzml).LoadAllStaticData();
            var ms1 = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1).OrderBy(s => s.OneBasedScanNumber).ToList();
            var slice = ms1.Skip(Math.Max(0, ms1.Count / 2 - 75)).Take(150).ToList();

            var perScan = new List<IReadOnlyList<IsotopicEnvelope>>();
            foreach (var s in slice) perScan.Add(Deconvoluter.Deconvolute(s.MassSpectrum, p).ToList());

            var tracer = new MassSpectrometry.Deconvolution.FeatureTracing.MetaFlashDeconMassFeatureTracer();
            var (input, traces) = tracer.DetectTracesDiagnostic(slice, perScan);

            var sb = new System.Text.StringBuilder();
            foreach (var (scanIdx, rtSec, peaks) in input)
            {
                sb.Append($"SPEC\t{scanIdx}\t{rtSec.ToString("R", CultureInfo.InvariantCulture)}\t{peaks.Count}\n");
                foreach (var (mass, intensity) in peaks)
                    sb.Append(mass.ToString("R", CultureInfo.InvariantCulture)).Append('\t').Append(intensity.ToString("R", CultureInfo.InvariantCulture)).Append('\n');
            }
            File.WriteAllText(MtdPeaks, sb.ToString());

            var sbo = new System.Text.StringBuilder();
            foreach (var (size, rtMin, rtMax, centroid) in traces)
                sbo.Append($"TRACE\t{size}\t{rtMin.ToString("R", CultureInfo.InvariantCulture)}\t{rtMax.ToString("R", CultureInfo.InvariantCulture)}\t{centroid.ToString("R", CultureInfo.InvariantCulture)}\n");
            File.WriteAllText(MtdOurs, sbo.ToString());
            TestContext.Progress.WriteLine($"ours traces={traces.Count}, scans={input.Count} -> {MtdPeaks}. Run masstrace_snip, then MassTrace_VsOpenMS.");
        }

        [Test]
        public void MassTrace_VsOpenMS()
        {
            Assume.That(File.Exists(MtdOurs), $"missing {MtdOurs} (run MassTrace_DumpForOpenMS)");
            Assume.That(File.Exists(MtdCpp), $"missing {MtdCpp} (run masstrace_snip)");
            (int size, double rtMin, double rtMax, double centroid) Parse(string l)
            { var t = l.Split('\t'); return (int.Parse(t[1]), double.Parse(t[2], CultureInfo.InvariantCulture), double.Parse(t[3], CultureInfo.InvariantCulture), double.Parse(t[4], CultureInfo.InvariantCulture)); }
            var ours = File.ReadLines(MtdOurs).Where(l => l.StartsWith("TRACE")).Select(Parse).OrderBy(x => x.centroid).ToList();
            var cpp = File.ReadLines(MtdCpp).Where(l => l.StartsWith("TRACE")).Select(Parse).OrderBy(x => x.centroid).ToList();

            // a trace matches if centroid within 10 ppm, RT ranges overlap, and size is identical
            bool Match((int size, double rtMin, double rtMax, double centroid) a, (int size, double rtMin, double rtMax, double centroid) b)
                => Math.Abs(a.centroid - b.centroid) <= b.centroid * 10e-6
                   && a.rtMin <= b.rtMax && b.rtMin <= a.rtMax && a.size == b.size;
            int MatchCount(List<(int,double,double,double)> A, List<(int,double,double,double)> B)
            { int n = 0; foreach (var a in A) if (B.Any(b => Match(a, b))) n++; return n; }

            int oursMatched = MatchCount(ours, cpp);
            int cppMatched = MatchCount(cpp, ours);
            TestContext.Progress.WriteLine($"ours traces={ours.Count}  openms traces={cpp.Count}");
            TestContext.Progress.WriteLine($"ours matching an OpenMS trace (centroid 10ppm + RT overlap + same size): {oursMatched}/{ours.Count}");
            TestContext.Progress.WriteLine($"OpenMS traces matched by ours: {cppMatched}/{cpp.Count}");
        }

        [Test]
        public void Collapse_MergeStats()
        {
            Assume.That(File.Exists(Mzml), $"missing {Mzml}");
            var p = new MetaFlashDeconParameters(minCharge: 1, maxCharge: 60);
            var dataFile = MsDataFileReader.GetDataFile(Mzml).LoadAllStaticData();
            var ms1 = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1).OrderBy(s => s.OneBasedScanNumber).ToList();
            var slice = ms1.Skip(Math.Max(0, ms1.Count / 2 - 30)).Take(60).ToList();
            var perScan = new List<IReadOnlyList<IsotopicEnvelope>>();
            foreach (var s in slice) perScan.Add(Deconvoluter.Deconvolute(s.MassSpectrum, p).ToList());
            var tracer = new MassSpectrometry.Deconvolution.FeatureTracing.MetaFlashDeconMassFeatureTracer();
            var (input, _) = tracer.DetectTracesDiagnostic(slice, perScan);
            int env = perScan.Sum(e => e.Count);
            int peaks = input.Sum(i => i.peaks.Count);
            TestContext.Progress.WriteLine($"COLLAPSE: scans={slice.Count} envelopes={env} collapsedPeaks={peaks} mergedAway={env - peaks} ({(env>0?(double)(env-peaks)/env:0):P2})");
        }

        // exactOnly: match within ppm; else allow ±1 isotope (Da) offset within ppm of the shifted mass.
        private static int FindMatch(List<double> real, double mass, double ppm, bool[] used, bool exactOnly)
        {
            int best = -1; double bestErr = double.MaxValue;
            for (int i = 0; i < real.Count; i++)
            {
                if (used[i]) continue;
                double tolDa = real[i] * ppm * 1e-6;
                double d = Math.Abs(real[i] - mass);
                if (!exactOnly) d = Math.Min(d, Math.Min(Math.Abs(real[i] - (mass + IsoDa)), Math.Abs(real[i] - (mass - IsoDa))));
                double allow = exactOnly ? tolDa : tolDa + 0.02;
                if (d <= allow && d < bestErr) { bestErr = d; best = i; }
            }
            return best;
        }

        private static Dictionary<int, List<double>> ParseMsalign(string path, out Dictionary<int, double> rtSeconds)
        {
            var result = new Dictionary<int, List<double>>();
            rtSeconds = new Dictionary<int, double>();
            int scan = -1; List<double> masses = null; double rt = -1;
            foreach (var raw in File.ReadLines(path))
            {
                var line = raw.Trim();
                if (line == "BEGIN IONS") { scan = -1; masses = new List<double>(); rt = -1; }
                else if (line.StartsWith("SCANS=")) scan = int.Parse(line.Substring(6), CultureInfo.InvariantCulture);
                else if (line.StartsWith("RETENTION_TIME=")) rt = double.Parse(line.Substring(15), CultureInfo.InvariantCulture);
                else if (line == "END IONS") { if (scan >= 0 && masses != null) { result[scan] = masses; rtSeconds[scan] = rt; } }
                else if (masses != null && line.Length > 0 && char.IsDigit(line[0]))
                {
                    var t = line.Split(new[] { '\t', ' ' }, StringSplitOptions.RemoveEmptyEntries);
                    if (t.Length >= 1 && double.TryParse(t[0], NumberStyles.Any, CultureInfo.InvariantCulture, out double m))
                        masses.Add(m);
                }
            }
            return result;
        }
    }
}
