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
