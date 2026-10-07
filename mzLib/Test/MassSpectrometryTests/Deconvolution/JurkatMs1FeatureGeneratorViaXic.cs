using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Deconvolution.Consensus;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.MassSpectrometryTests.Deconvolution
{
    /// <summary>
    /// YIELD-COMPARISON variant of <see cref="JurkatMs1FeatureGenerator"/>: produces the same
    /// <c>_ms1.feature</c> output, but does the per-scan grouping with mzLib's EXISTING XIC engine
    /// (<see cref="MassIndexingEngine"/> + <see cref="IndexingEngine{T}.GetAllXics"/>) instead of
    /// the PR's <see cref="MassTraceBuilder"/>. Everything downstream -- <see cref="TraceCorrector"/>
    /// (off-by-one + deamidation split) and <see cref="MassFeatureBuilder"/> (cross-charge stitch)
    /// and the writer -- is the IDENTICAL production code path. The per-envelope mass/intensity/charge
    /// are read straight off each XIC peak's wrapped <see cref="IndexedMass.IsotopicEnvelope"/>, so the
    /// only thing that differs from the consensus generator is which envelopes get grouped together.
    ///
    /// This answers nic's reuse question at the decisive metric: feed these features through the same
    /// MetaMorpheus FromFile-only search and compare top-down proteoforms@q&lt;0.01 to the 1,100 baseline.
    ///
    /// Output goes to a dedicated folder so it never collides with the consensus generator's output.
    /// Run: dotnet test --filter FullyQualifiedName~JurkatMs1FeatureGeneratorViaXic
    /// </summary>
    [TestFixture]
    [Explicit("Generates XIC-grouped Ms1Feature files for the local Jurkat TopDown dataset; needs E:\\Projects\\JurkatTopDown")]
    [ExcludeFromCodeCoverage]
    public class JurkatMs1FeatureGeneratorViaXic
    {
        private const string FolderPath = @"E:\Projects\JurkatTopDown";
        private const string OutFolder = @"E:\Projects\JurkatTopDown\XicVariant";

        // Top-down Classic-decon parameters (identical to the consensus generator).
        private const int MinCharge = 1;
        private const int MaxCharge = 60;
        private const double DeconPpm = 4.0;
        private const double IntensityRatio = 3.0;

        // Tracing parameters (match the consensus generator's window so only the ALGORITHM differs).
        private const double TraceToleranceDa = 1.5;
        private const int MaxGap = 1;
        private const double CrossChargeMassPpm = 10.0;

        // Above the consensus 30 kDa default so the bin cap is not the differentiator.
        private const int MaxMassDa = 60000;

        [Test]
        public void GenerateXicGroupedMs1FeatureFilesForAllJurkatMzmls()
        {
            Assume.That(Directory.Exists(FolderPath), $"Folder not found: {FolderPath}");
            Directory.CreateDirectory(OutFolder);
            var mzmls = Directory.EnumerateFiles(FolderPath, "*.mzML").OrderBy(p => p).ToList();
            Assume.That(mzmls, Is.Not.Empty, $"No mzML files at {FolderPath}");

            TestContext.Progress.WriteLine($"[XIC] Found {mzmls.Count} mzMLs in {FolderPath}");
            var overall = Stopwatch.StartNew();
            int produced = 0, skipped = 0, failed = 0;
            var failures = new List<string>();

            foreach (string mzml in mzmls)
            {
                string name = Path.GetFileNameWithoutExtension(mzml);
                string outputPath = Path.Combine(OutFolder, name + "_ms1.feature");
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
                    var scanArray = ms1Scans.ToArray();

                    var deconParams = new ClassicDeconvolutionParameters(
                        minCharge: MinCharge, maxCharge: MaxCharge,
                        deconPpm: DeconPpm, intensityRatio: IntensityRatio);

                    // EXISTING engine: deconvolute + index once, then trace all XICs (charge-locked,
                    // gap-tolerant, same absolute-Da window as the consensus tracer).
                    var engine = new MassIndexingEngine { MaxMass = MaxMassDa };
                    if (!engine.IndexPeaks(scanArray, deconParams))
                        throw new MzLibException("IndexPeaks returned false (no indexable scans).");
                    var xics = engine.GetAllXics(
                        new AbsoluteTolerance(TraceToleranceDa), MaxGap, double.MaxValue, numPeakThreshold: 1);

                    // Adapt each XIC into the MassTrace the downstream engines expect, reading the
                    // exact envelope mass/intensity/charge off the wrapped IsotopicEnvelope.
                    var traces = new List<MassTrace>(xics.Count);
                    int id = 1;
                    foreach (var xic in xics)
                        traces.Add(XicToMassTrace(xic, id++, ms1Scans));

                    var corrected = traces.SelectMany(t => TraceCorrector.Correct(t)).ToList();
                    var features = MassFeatureBuilder.BuildFeatures(corrected, CrossChargeMassPpm);

                    Ms1FeatureFile.FromMassFeatures(features).WriteResults(outputPath);

                    sw.Stop();
                    produced++;
                    TestContext.Progress.WriteLine(
                        $"[{name}] done: xics={xics.Count} traces={traces.Count} features={features.Count} in {sw.Elapsed.TotalSeconds:F0}s");
                }
                catch (Exception ex)
                {
                    failed++;
                    failures.Add($"{name}: {ex.GetType().Name} {ex.Message}");
                    TestContext.Progress.WriteLine($"[FAIL] {name}: {ex.Message}");
                }
            }

            overall.Stop();
            TestContext.Progress.WriteLine(
                $"[XIC] produced={produced} skipped={skipped} failed={failed} in {overall.Elapsed.TotalMinutes:F1} min");
            if (failures.Count > 0)
                TestContext.Progress.WriteLine("Failures:\n  " + string.Join("\n  ", failures));
            Assert.That(failed, Is.EqualTo(0), "One or more files failed to generate; see log.");
        }

        /// <summary>
        /// Build the downstream <see cref="MassTrace"/> from an XIC. Charge/mass/intensity come from
        /// each peak's wrapped <see cref="IndexedMass.IsotopicEnvelope"/>, so the values are identical
        /// to what the consensus generator feeds; only the grouping (which envelopes are together) differs.
        /// </summary>
        private static MassTrace XicToMassTrace(ExtractedIonChromatogram xic, int id, List<MsDataScan> ms1Scans)
        {
            var peaks = xic.Peaks.OfType<IndexedMass>().OrderBy(p => p.ZeroBasedScanIndex).ToList();
            var mt = new MassTrace
            {
                Id = id,
                Charge = peaks.Count > 0 ? peaks[0].IsotopicEnvelope.Charge : 0,
                AnchorMass = peaks.Count > 0 ? peaks[0].IsotopicEnvelope.MonoisotopicMass : 0,
            };
            foreach (var p in peaks)
            {
                var env = p.IsotopicEnvelope;
                mt.Envelopes.Add((
                    p.ZeroBasedScanIndex,
                    ms1Scans[p.ZeroBasedScanIndex].OneBasedScanNumber,
                    p.RetentionTime,
                    env.MonoisotopicMass,
                    env.TotalIntensity));
            }
            return mt;
        }
    }
}
