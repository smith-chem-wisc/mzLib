// TestFLASHDeconvMzmlToMsAlign.cs
//
// PURPOSE: Integration test that loads a real mzML file from disk, runs FLASHDeconvolution
// on every MS1 scan, and writes the results to a _ms1.msalign file inside a subfolder
// named after the input file (without extension).
//
// HOW TO USE:
//   1. Set MzmlPath to the full path of your mzML file (or use the [TestCase] attribute).
//   2. Run the test. The output subfolder and msalign file are created next to the mzML.
//   3. Compare the output with the reference FLASHDeconv msalign.
//
// EXAMPLE:
//   Input:  E:\Projects\...\190226_FIlg_2_FD_500ng.mzML
//   Output: E:\Projects\...\190226_FIlg_2_FD_500ng\190226_FIlg_2_FD_500ng_ms1.msalign

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
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
        // Change this path (or add [TestCase(...)] entries) to point at your mzML.
        private const string MzmlPath =
            @"E:\Projects\Deconvolution\FlashDeconv_MSV000084001\updates\2019-12-04_kyowonjeong_e9e915a8\raw\190226_FIlg_2_FD_500ng.mzML";

        // ────────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Loads the mzML at <see cref="MzmlPath"/>, deconvolutes every MS1 scan
        /// using FLASHDeconvolution, and writes a _ms1.msalign file into a subfolder
        /// that is named after the mzML file (without its extension).
        /// </summary>
        [Test]
        [TestCase(MzmlPath)]
        public static void DeconvoluteMzmlAndWriteMsAlign(string mzmlPath)
        {
            // ── Pre-condition checks ─────────────────────────────────────────────
            Assert.That(File.Exists(mzmlPath),
                $"mzML file not found: {mzmlPath}\n" +
                "Set MzmlPath at the top of TestFLASHDeconvMzmlToMsAlign.cs to a valid path.");

            // ── Build output paths ───────────────────────────────────────────────
            string mzmlDir = Path.GetDirectoryName(mzmlPath)!;
            string baseName = Path.GetFileNameWithoutExtension(mzmlPath);   // e.g. 190226_FIlg_2_FD_500ng
            string subfolderPath = Path.Combine(mzmlDir, baseName);
            string outputPath = Path.Combine(subfolderPath, $"mzLib_{baseName}_ms1.msalign");

            Directory.CreateDirectory(subfolderPath);
            Console.WriteLine($"[Output folder] {subfolderPath}");
            Console.WriteLine($"[Output file  ] {outputPath}");

            // ── Load mzML ────────────────────────────────────────────────────────
            Console.WriteLine($"Loading {mzmlPath} ...");
            MsDataFile dataFile = MsDataFileReader.GetDataFile(mzmlPath).LoadAllStaticData();
            List<MsDataScan> ms1Scans = dataFile.GetAllScansList()
                                                .Where(s => s.MsnOrder == 1)
                                                .ToList();
            Console.WriteLine($"  Total scans loaded : {dataFile.NumSpectra}");
            Console.WriteLine($"  MS1 scans          : {ms1Scans.Count}");
            Assert.That(ms1Scans.Count, Is.GreaterThan(0), "No MS1 scans found in the mzML file.");

            // ── Deconvolute ──────────────────────────────────────────────────────
            var deconParams = new FLASHDeconvolutionParameters();

            Console.WriteLine("Deconvoluting MS1 scans ...");
            var results = new List<(MsDataScan Scan, List<IsotopicEnvelope> Envelopes)>();
            int scansWithResults = 0;
            int totalEnvelopes = 0;

            foreach (MsDataScan scan in ms1Scans)
            {
                List<IsotopicEnvelope> envelopes = Deconvoluter
                    .Deconvolute(scan, deconParams)
                    .ToList();

                results.Add((scan, envelopes));

                if (envelopes.Count > 0)
                {
                    scansWithResults++;
                    totalEnvelopes += envelopes.Count;
                }
            }

            Console.WriteLine($"  MS1 scans with ≥1 envelope : {scansWithResults}");
            Console.WriteLine($"  Total isotopic envelopes   : {totalEnvelopes}");

            // ── Write msalign ────────────────────────────────────────────────────
            Console.WriteLine($"Writing msalign to {outputPath} ...");
            WriteMsAlign(results, outputPath);

            // ── Verify output ────────────────────────────────────────────────────
            Assert.That(File.Exists(outputPath), $"msalign file was not created at {outputPath}");

            long fileSize = new FileInfo(outputPath).Length;
            Console.WriteLine($"  msalign file size : {fileSize:N0} bytes");
            Assert.That(fileSize, Is.GreaterThan(0), "msalign file was created but is empty.");

            // Quick sanity: count BEGIN IONS lines
            int beginIonsCount = File.ReadLines(outputPath)
                                     .Count(l => l.TrimStart().StartsWith("BEGIN IONS"));
            Console.WriteLine($"  BEGIN IONS blocks  : {beginIonsCount}");
            Assert.That(beginIonsCount, Is.EqualTo(scansWithResults),
                "Number of BEGIN IONS blocks does not match number of scans with results.");

            Console.WriteLine("Done. msalign written successfully.");
            Console.WriteLine($"\nNext step: compare\n  {outputPath}\nwith the reference FLASHDeconv output.");
        }

        // ── msAlign writer ───────────────────────────────────────────────────────
        // Inline here so the test file is self-contained.
        // Format follows: NewDeconvolutionMethod_ImplementationPrompt.md §4
        // and the existing MsAlign reader in Readers/ExternalFileTypes/MsAlign.cs.

        private static void WriteMsAlign(
            IEnumerable<(MsDataScan Scan, List<IsotopicEnvelope> Envelopes)> results,
            string outputPath)
        {
            using var writer = new StreamWriter(outputPath);

            foreach (var (scan, envelopes) in results)
            {
                if (envelopes == null || envelopes.Count == 0)
                    continue;

                double retentionTimeSec = scan.RetentionTime * 60.0;

                writer.WriteLine("BEGIN IONS");
                writer.WriteLine($"ID={scan.OneBasedScanNumber}");
                writer.WriteLine($"FRACTION_ID=0");
                writer.WriteLine($"SCANS={scan.OneBasedScanNumber}");
                writer.WriteLine($"RETENTION_TIME={retentionTimeSec:F2}");
                writer.WriteLine($"LEVEL=1");

                foreach (IsotopicEnvelope env in envelopes.OrderBy(e => e.MonoisotopicMass))
                {
                    int absCharge = Math.Abs(env.Charge);
                    writer.WriteLine($"{env.MonoisotopicMass:F6}\t{env.TotalIntensity:F2}\t{absCharge}");
                }

                writer.WriteLine("END IONS");
                writer.WriteLine();
            }
        }
    }
}