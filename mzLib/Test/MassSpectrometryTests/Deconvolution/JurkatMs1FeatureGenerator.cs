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
    /// EXPLICIT, throwaway driver: for every .mzML in <see cref="FolderPath"/>, runs the full
    /// consensus pipeline (Classic decon per MS1 -> trace grouping -> resolution-aware off-by-one
    /// correction with deamidation split -> cross-charge stitching) and emits a sibling
    /// <c>_ms1.feature</c> file via <see cref="Ms1FeatureFile.FromMassFeatures"/>.
    ///
    /// Top-down parameters (charge 1..60). Skip-if-exists so resuming after a failure is safe.
    /// Not part of CI; run manually with:
    ///   dotnet test --filter FullyQualifiedName~JurkatMs1FeatureGenerator
    /// </summary>
    [TestFixture]
    [Explicit("Generates Ms1Feature files for the local Jurkat TopDown dataset; needs E:\\Projects\\JurkatTopDown")]
    [ExcludeFromCodeCoverage]
    public class JurkatMs1FeatureGenerator
    {
        private const string FolderPath = @"E:\Projects\JurkatTopDown";

        // Top-down Classic-decon parameters.
        private const int MinCharge = 1;
        private const int MaxCharge = 60;
        private const double DeconPpm = 4.0;
        private const double IntensityRatio = 3.0;

        // Consensus-tracing parameters (match the end-to-end driver).
        private const double PassBToleranceDa = 1.5;
        private const int MaxGap = 1;
        private const double CrossChargeMassPpm = 10.0;

        [Test]
        public void GenerateMs1FeatureFilesForAllJurkatMzmls()
        {
            Assume.That(Directory.Exists(FolderPath), $"Folder not found: {FolderPath}");
            var mzmls = Directory.EnumerateFiles(FolderPath, "*.mzML").OrderBy(p => p).ToList();
            Assume.That(mzmls, Is.Not.Empty, $"No mzML files at {FolderPath}");

            TestContext.Progress.WriteLine($"[Jurkat] Found {mzmls.Count} mzMLs in {FolderPath}");
            var overall = Stopwatch.StartNew();

            int produced = 0, skipped = 0, failed = 0;
            var failures = new List<string>();

            foreach (string mzml in mzmls)
            {
                string name = Path.GetFileNameWithoutExtension(mzml);
                string outputPath = Path.Combine(FolderPath, name + "_ms1.feature");
                if (File.Exists(outputPath))
                {
                    TestContext.Progress.WriteLine($"[skip] {name}: output already exists");
                    skipped++;
                    continue;
                }

                try
                {
                    var sw = Stopwatch.StartNew();
                    TestContext.Progress.WriteLine($"[{name}] starting...");

                    var dataFile = MsDataFileReader.GetDataFile(mzml).LoadAllStaticData();
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
                        PassBToleranceDa,
                        MaxGap);
                    // Correct() now returns List<CorrectedTrace> (a trace may split into a distinct
                    // co-grouped species, e.g. deamidation), so flatten with SelectMany.
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
            TestContext.Progress.WriteLine(
                $"[Jurkat] done: produced={produced}, skipped={skipped}, failed={failed}, total elapsed {overall.Elapsed}");
            if (failures.Count > 0)
            {
                foreach (var f in failures) TestContext.Progress.WriteLine("  fail: " + f);
                Assert.Fail($"{failures.Count} mzML(s) failed; see output.");
            }
        }
    }
}
