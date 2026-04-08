// TestFLASHDeconvMzmlToMsAlign.cs
//
// PURPOSE: Integration test that loads a real mzML file from disk, runs FLASHDeconvolution
// on every MS1 scan, and writes:
//   1. A _ms1.msalign file (masses, intensities, charges per scan)
//   2. A _ms1.tsv file    (per-envelope: Qscore, SimpleQvalue, ppm error, etc.)
// Both files are written into a subfolder named after the input mzML file.
//
// PIPELINE:
//   Deconvoluter.Deconvolute()  — internally runs Steps 1-5, dedup, Qscoring
//                                 returns envelopes with Score = Qscore in [0,1]
//
// NOTE: Score on returned envelopes is the Qscore, not the original cosine.
// The IsotopeCosine column in the TSV is therefore unavailable without modifying
// the return type of Deconvolute(). The column is written as -1 as a placeholder.
// The Qscore is computed inside the algorithm using the exact ppm error from
// Step 3 (isotope index known during recruitment), which is more accurate than
// the post-hoc approximation.
//
// NOTE ON SimpleQvalue:
//   OpenMS computes q-values via two extra decoy deconvolution passes per spectrum.
//   We do not implement that here. SimpleQvalue = 1 - Qscore, monotone enforced —
//   an honest approximation using the logistic regression Qscore as a local FDR.

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.Deconvolution
{
    [TestFixture]
    public class TestFLASHDeconvMzmlToMsAlign
    {
        // ── Configure your input file here ──────────────────────────────────────
        private const string MzmlPath =
            @"E:\Projects\Deconvolution\FlashDeconv_MSV000084001\updates\2019-12-04_kyowonjeong_e9e915a8\raw\190226_FIlg_2_FD_500ng.mzML";
        // ────────────────────────────────────────────────────────────────────────

        [Test]
        [TestCase(MzmlPath)]
        public static void DeconvoluteMzmlAndWriteOutputs(string mzmlPath)
        {
            // ── Pre-condition ────────────────────────────────────────────────────
            Assert.That(File.Exists(mzmlPath), $"mzML file not found: {mzmlPath}");

            // ── Output paths ─────────────────────────────────────────────────────
            string mzmlDir = Path.GetDirectoryName(mzmlPath)!;
            string baseName = Path.GetFileNameWithoutExtension(mzmlPath);
            string subfolderPath = Path.Combine(mzmlDir, baseName);
            string msalignPath = Path.Combine(subfolderPath, $"mzLib_{baseName}_ms1.msalign");
            string tsvPath = Path.Combine(subfolderPath, $"mzLib_{baseName}_ms1.tsv");

            Directory.CreateDirectory(subfolderPath);
            Console.WriteLine($"[Output folder] {subfolderPath}");

            // ── Load mzML ────────────────────────────────────────────────────────
            Console.WriteLine($"Loading {mzmlPath} ...");
            MsDataFile dataFile = MsDataFileReader.GetDataFile(mzmlPath).LoadAllStaticData();
            List<MsDataScan> ms1Scans = dataFile.GetAllScansList()
                                                .Where(s => s.MsnOrder == 1)
                                                .ToList();
            Console.WriteLine($"  Total scans : {dataFile.NumSpectra}");
            Console.WriteLine($"  MS1 scans   : {ms1Scans.Count}");
            Assert.That(ms1Scans.Count, Is.GreaterThan(0), "No MS1 scans found.");

            // ── Deconvolute ──────────────────────────────────────────────────────
            // Deconvolute() runs the full pipeline internally:
            //   Steps 1-5 → deduplication → Qscoring
            // Returned envelopes have Score = Qscore in [0,1].
            var deconParams = new FLASHDeconvolutionParameters();
            Console.WriteLine("Deconvoluting ...");

            var msalignResults = new List<(MsDataScan Scan, List<IsotopicEnvelope> Envelopes)>();
            var tsvResults = new List<(MsDataScan Scan, IsotopicEnvelope Envelope)>();

            int scansWithResults = 0;
            int totalEnvelopes = 0;

            foreach (MsDataScan scan in ms1Scans)
            {
                // Score on returned envelopes = Qscore (set internally by FLASHDeconvScorer)
                List<IsotopicEnvelope> envelopes = Deconvoluter
                    .Deconvolute(scan, deconParams)
                    .ToList();

                msalignResults.Add((scan, envelopes));

                foreach (var env in envelopes)
                    tsvResults.Add((scan, env));

                if (envelopes.Count > 0)
                {
                    scansWithResults++;
                    totalEnvelopes += envelopes.Count;
                }
            }

            Console.WriteLine($"  Scans with ≥1 envelope : {scansWithResults}");
            Console.WriteLine($"  Total envelopes        : {totalEnvelopes}");

            // ── Write msalign ────────────────────────────────────────────────────
            Console.WriteLine($"Writing {msalignPath} ...");
            WriteMsAlign(msalignResults, msalignPath);
            Assert.That(File.Exists(msalignPath));
            Assert.That(new FileInfo(msalignPath).Length, Is.GreaterThan(0));
            Console.WriteLine($"  size : {new FileInfo(msalignPath).Length:N0} bytes");

            // ── Write TSV ────────────────────────────────────────────────────────
            Console.WriteLine($"Writing {tsvPath} ...");
            WriteTsv(tsvResults, tsvPath);
            Assert.That(File.Exists(tsvPath));
            Assert.That(new FileInfo(tsvPath).Length, Is.GreaterThan(0));
            Console.WriteLine($"  size : {new FileInfo(tsvPath).Length:N0} bytes");

            Console.WriteLine("Done.");
        }

        // ════════════════════════════════════════════════════════════════════════
        // msalign writer
        // ════════════════════════════════════════════════════════════════════════

        private static void WriteMsAlign(
            IEnumerable<(MsDataScan Scan, List<IsotopicEnvelope> Envelopes)> results,
            string outputPath)
        {
            using var writer = new StreamWriter(outputPath);
            foreach (var (scan, envelopes) in results)
            {
                if (envelopes == null || envelopes.Count == 0) continue;
                writer.WriteLine("BEGIN IONS");
                writer.WriteLine($"ID={scan.OneBasedScanNumber}");
                writer.WriteLine($"FRACTION_ID=0");
                writer.WriteLine($"SCANS={scan.OneBasedScanNumber}");
                writer.WriteLine($"RETENTION_TIME={scan.RetentionTime * 60.0:F2}");
                writer.WriteLine($"LEVEL=1");
                foreach (IsotopicEnvelope env in envelopes.OrderBy(e => e.MonoisotopicMass))
                    writer.WriteLine($"{env.MonoisotopicMass:F6}\t{env.TotalIntensity:F2}\t{Math.Abs(env.Charge)}");
                writer.WriteLine("END IONS");
                writer.WriteLine();
            }
        }

        // ════════════════════════════════════════════════════════════════════════
        // TSV writer
        // Columns mirror OpenMS FLASHDeconv _ms1.tsv where possible.
        // IsotopeCosine = -1 (not available; Score has been replaced by Qscore).
        // SimpleQvalue  = 1 - Qscore, monotone non-decreasing as Qscore decreases.
        // ════════════════════════════════════════════════════════════════════════

        private static readonly string[] TsvHeader =
        {
            "ScanNum", "RetentionTime", "MonoisotopicMass", "AverageMass",
            "SumIntensity", "MinCharge", "MaxCharge", "PeakCount",
            "IsotopeCosine", "Qscore", "SimpleQvalue", "RepresentativeCharge", "AvgPPMError"
        };

        private static void WriteTsv(
            List<(MsDataScan Scan, IsotopicEnvelope Envelope)> rows,
            string outputPath)
        {
            if (rows.Count == 0) return;

            // SimpleQvalue = 1 − Qscore, monotone enforced globally
            var simpleQvalue = new double[rows.Count];
            double runningMin = 1.0;
            foreach (var (r, i) in rows.Select((r, i) => (r, i))
                                       .OrderByDescending(x => x.r.Envelope.Score))
            {
                runningMin = Math.Min(runningMin, 1.0 - r.Envelope.Score);
                simpleQvalue[i] = runningMin;
            }
            double runningMax = 0.0;
            foreach (var (r, i) in rows.Select((r, i) => (r, i))
                                       .OrderBy(x => x.r.Envelope.Score))
            {
                simpleQvalue[i] = Math.Max(simpleQvalue[i], runningMax);
                runningMax = simpleQvalue[i];
            }

            using var writer = new StreamWriter(outputPath);
            writer.WriteLine(string.Join("\t", TsvHeader));

            foreach (var (r, i) in rows.Select((r, i) => (r, i))
                                       .OrderBy(x => x.r.Scan.OneBasedScanNumber)
                                       .ThenByDescending(x => x.r.Envelope.Score))
            {
                var (scan, env) = r;
                int absCharge = Math.Abs(env.Charge);
                double apexMz = env.Peaks.Count > 0
                    ? env.Peaks.MaxBy(p => p.intensity).mz
                    : env.MonoisotopicMass.ToMz(env.Charge);
                double averageMass = apexMz * absCharge;
                double avgPpmError = RecomputeAvgPpmError(env);

                writer.WriteLine(string.Join("\t", new object[]
                {
                    scan.OneBasedScanNumber,
                    $"{scan.RetentionTime * 60.0:F2}",
                    $"{env.MonoisotopicMass:F6}",
                    $"{averageMass:F6}",
                    $"{env.TotalIntensity:F2}",
                    absCharge, absCharge,
                    env.Peaks.Count,
                    "-1",                          // IsotopeCosine: not available (Score = Qscore)
                    $"{env.Score:F6}",             // Qscore
                    $"{simpleQvalue[i]:F6}",
                    absCharge,
                    $"{avgPpmError:F4}"
                }));
            }
        }

        // Post-hoc ppm error for display in TSV (approximate — the exact value
        // was used inside Deconvolute() for Qscore computation but is not returned).
        private static double RecomputeAvgPpmError(IsotopicEnvelope env)
        {
            if (env.Peaks.Count == 0 || env.Charge == 0) return 0.0;
            int absCharge = Math.Abs(env.Charge);
            double isotopeStepMz = Constants.C13MinusC12 / absCharge;
            double apexMz = env.Peaks.MaxBy(p => p.intensity).mz;
            double totalPpm = 0.0;
            foreach (var (mz, _) in env.Peaks)
            {
                int n = (int)Math.Round((mz - apexMz) / isotopeStepMz);
                double theorMz = apexMz + n * isotopeStepMz;
                totalPpm += Math.Abs(mz - theorMz) / theorMz * 1e6;
            }
            return totalPpm / env.Peaks.Count;
        }
    }
}