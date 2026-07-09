using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Deconvolution.Consensus;
using NUnit.Framework;
using Readers;

namespace Test.MassSpectrometryTests.Deconvolution
{
    /// <summary>
    /// EXPLICIT, throwaway driver: like <see cref="JurkatMs1FeatureGenerator"/> but for the
    /// bottom-up HEK293 Velos raws. Globs *.raw (the 4 full files; the _snip_*.mzML are not
    /// matched) and emits a sibling <c>_ms1.feature</c> per raw via the new deamidation-split
    /// pipeline. Bottom-up charge range (1..12). Skip-if-exists; not part of CI.
    ///   dotnet test --filter FullyQualifiedName~VelosBuMs1FeatureGenerator
    /// </summary>
    [TestFixture]
    [Explicit("Generates Ms1Feature files for the local HEK293 Velos BU raws; needs E:\\TestData\\MetaMorpheus\\BottomUp")]
    [ExcludeFromCodeCoverage]
    public class VelosBuMs1FeatureGenerator
    {
        private const string FolderPath = @"E:\TestData\MetaMorpheus\BottomUp";

        // Bottom-up Classic-decon parameters (tryptic peptides rarely exceed +6; 12 matches
        // MetaMorpheus's BU DeconvolutionMaxAssumedChargeState default).
        private const int MinCharge = 1;
        private const int MaxCharge = 12;
        private const double DeconPpm = 4.0;
        private const double IntensityRatio = 3.0;

        private const double PassBToleranceDa = 1.5;
        private const int MaxGap = 1;
        private const double CrossChargeMassPpm = 10.0;

        [Test]
        public void GenerateMs1FeatureFilesForAllVelosRaws()
        {
            Assume.That(Directory.Exists(FolderPath), $"Folder not found: {FolderPath}");
            var raws = Directory.EnumerateFiles(FolderPath, "*.raw").OrderBy(p => p).ToList();
            Assume.That(raws, Is.Not.Empty, $"No .raw files at {FolderPath}");

            TestContext.Progress.WriteLine($"[Velos BU] Found {raws.Count} raws in {FolderPath}");
            var overall = Stopwatch.StartNew();
            int produced = 0, skipped = 0, failed = 0;
            var failures = new List<string>();

            foreach (string raw in raws)
            {
                string name = Path.GetFileNameWithoutExtension(raw);
                string outputPath = Path.Combine(FolderPath, name + "_ms1.feature");
                if (File.Exists(outputPath)) { TestContext.Progress.WriteLine($"[skip] {name}"); skipped++; continue; }

                try
                {
                    var sw = Stopwatch.StartNew();
                    TestContext.Progress.WriteLine($"[{name}] starting...");

                    var dataFile = MsDataFileReader.GetDataFile(raw).LoadAllStaticData();
                    var ms1Scans = dataFile.GetAllScansList()
                        .Where(s => s.MsnOrder == 1)
                        .OrderBy(s => s.OneBasedScanNumber)
                        .ToList();

                    var deconParams = new ClassicDeconvolutionParameters(
                        minCharge: MinCharge, maxCharge: MaxCharge,
                        deconPpm: DeconPpm, intensityRatio: IntensityRatio);
                    var perScanEnvelopes = ms1Scans
                        .Select(s => Deconvoluter.Deconvolute(s.MassSpectrum, deconParams).ToList())
                        .ToList();
                    int totalEnv = perScanEnvelopes.Sum(s => s.Count);

                    var traces = MassTraceBuilder.BuildTraces(
                        ms1Scans,
                        perScanEnvelopes.Select(scan => (IReadOnlyList<IsotopicEnvelope>)scan).ToList(),
                        PassBToleranceDa, MaxGap);
                    var corrected = traces.SelectMany(t => TraceCorrector.Correct(t)).ToList();
                    var features = MassFeatureBuilder.BuildFeatures(corrected, CrossChargeMassPpm);

                    Ms1FeatureFile.FromMassFeatures(features).WriteResults(outputPath);
                    ConsensusFeatureScoring.WriteScores(outputPath + ".scores.tsv", features, perScanEnvelopes);

                    sw.Stop();
                    long bytes = new FileInfo(outputPath).Length;
                    TestContext.Progress.WriteLine(
                        $"[{name}]   MS1={ms1Scans.Count} env={totalEnv} traces={traces.Count} " +
                        $"features={features.Count} (multi-z={features.Count(f => f.ChargeCount >= 2)}) " +
                        $"wrote {bytes / 1024.0:F1} KB in {sw.Elapsed}");
                    produced++;
                }
                catch (Exception ex)
                {
                    failed++;
                    failures.Add($"{name}: {ex.GetType().Name}: {ex.Message}");
                    TestContext.Progress.WriteLine($"[{name}] FAILED: {ex.Message}");
                }
            }

            overall.Stop();
            TestContext.Progress.WriteLine($"[Velos BU] done: produced={produced}, skipped={skipped}, failed={failed}, elapsed {overall.Elapsed}");
            if (failures.Count > 0) { foreach (var f in failures) TestContext.Progress.WriteLine("  fail: " + f); Assert.Fail($"{failures.Count} raw(s) failed."); }
        }
    }
}
