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
        private const int SampleEvery = 80; // ~40 scans across the run

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
