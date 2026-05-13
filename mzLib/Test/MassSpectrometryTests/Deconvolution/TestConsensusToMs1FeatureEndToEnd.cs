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
    /// End-to-end driver: run the full consensus pipeline (Classic decon per
    /// scan -&gt; trace grouping -&gt; off-by-one correction -&gt; cross-charge feature
    /// stitching) against a local mzML and write the result as a FLASHDeconv-
    /// style <c>_ms1.feature</c> file via <see cref="Ms1FeatureFile.FromMassFeatures"/>.
    ///
    /// Lives as an <c>[Explicit]</c> test because it consumes a 165 MB local
    /// mzML at a hard-coded path that's not part of the repo and not part of CI.
    /// Run with:
    /// <code>
    /// dotnet test --filter FullyQualifiedName~TestConsensusToMs1FeatureEndToEnd
    /// </code>
    ///
    /// Pipeline parameters match the values established in the consensus-
    /// tracing research arc on the yeast top-down data: Pass B grouping at
    /// 1.5 Da (intentionally loose to cohort off-by-one twins for correction),
    /// 1-scan gap tolerance, 10 ppm cross-charge mass agreement.
    /// </summary>
    [TestFixture]
    [Explicit("End-to-end consensus pipeline driver; needs local mzML at E:\\TestData\\MetaMorpheus\\.")]
    [ExcludeFromCodeCoverage]
    public class TestConsensusToMs1FeatureEndToEnd
    {
        // Local-only paths.
        private const string SourceMzMl =
            @"E:\TestData\MetaMorpheus\05-26-17_B7A_yeast_td_fract7_rep1.mzML";
        private const string OutputFeatureFile =
            @"E:\TestData\MetaMorpheus\05-26-17_B7A_yeast_td_fract7_rep1_ms1.feature";

        // Pipeline parameters (per NOTES.md big-snip defaults on the consensus-
        // tracing branch). Classic decon at top-down-appropriate charge range.
        private const double PassBToleranceDa = 1.5;
        private const int    MaxGap            = 1;
        private const double CrossChargeMassPpm = 10.0;

        [Test]
        public void RunFullPipelineAndEmitMs1FeatureFile()
        {
            Assume.That(File.Exists(SourceMzMl),
                $"mzML not found at {SourceMzMl}; install local test data or override the path.");

            var sw = Stopwatch.StartNew();

            // 1. Load the mzML and extract MS1 scans in order.
            var dataFile = MsDataFileReader.GetDataFile(SourceMzMl).LoadAllStaticData();
            var ms1Scans = dataFile.GetAllScansList()
                .Where(s => s.MsnOrder == 1)
                .OrderBy(s => s.OneBasedScanNumber)
                .ToList();
            Assume.That(ms1Scans, Is.Not.Empty);
            TestContext.Out.WriteLine($"MS1 scans loaded: {ms1Scans.Count}  (elapsed {sw.Elapsed})");

            // 2. Per-scan Classic decon. Top-down charge range (1..60) matches
            //    DeconvolutionMaxAssumedChargeState=60 in the MetaMorpheus task
            //    config used for the yeast TD benchmark.
            var deconParams = new ClassicDeconvolutionParameters(
                minCharge: 1,
                maxCharge: 60,
                deconPpm: 4.0,
                intensityRatio: 3.0);
            var perScanEnvelopes = ms1Scans
                .Select(s => Deconvoluter.Deconvolute(s.MassSpectrum, deconParams).ToList())
                .ToList();
            int totalEnvelopes = perScanEnvelopes.Sum(s => s.Count);
            TestContext.Out.WriteLine($"Classic envelopes total: {totalEnvelopes}  (elapsed {sw.Elapsed})");

            // 3. Pass B trace grouping (loose: 1.5 Da) so off-by-one twins
            //    cohort with their parents and the corrector can rescue them.
            var traces = MassTraceBuilder.BuildTraces(
                ms1Scans,
                perScanEnvelopes.Select(scan => (IReadOnlyList<IsotopicEnvelope>)scan).ToList(),
                PassBToleranceDa,
                MaxGap);
            TestContext.Out.WriteLine($"Traces built: {traces.Count}  (elapsed {sw.Elapsed})");

            // 4. Per-trace weighted-median off-by-one correction. Uniform
            //    weight matches NOTES.md's Phase 5 finding that weighting
            //    scheme barely changes consensus mass.
            var corrected = traces.Select(t => TraceCorrector.Correct(t)).ToList();
            int corrections = corrected.Sum(c => c.CorrectionCount);
            TestContext.Out.WriteLine($"Corrections applied: {corrections}  (elapsed {sw.Elapsed})");

            // 5. Cross-charge clustering.
            var features = MassFeatureBuilder.BuildFeatures(corrected, CrossChargeMassPpm);
            int multi = features.Count(f => f.ChargeCount >= 2);
            TestContext.Out.WriteLine($"MassFeatures: {features.Count} (multi-charge: {multi})  (elapsed {sw.Elapsed})");

            // 6. Write FLASHDeconv-style _ms1.feature via the new factory.
            Ms1FeatureFile.FromMassFeatures(features).WriteResults(OutputFeatureFile);

            sw.Stop();
            var bytes = new FileInfo(OutputFeatureFile).Length;
            TestContext.Out.WriteLine($"Wrote {bytes / 1024.0:F1} KB to {OutputFeatureFile}  (total elapsed {sw.Elapsed})");

            // Sanity: the file should round-trip cleanly.
            var roundTrip = FileReader.ReadFile<Ms1FeatureFile>(OutputFeatureFile);
            Assert.That(roundTrip.Results.Count, Is.EqualTo(features.Count),
                "Round-trip row count mismatch -- writer or reader regression.");
            TestContext.Out.WriteLine($"Round-trip read: {roundTrip.Results.Count} rows.");

            // Top 10 multi-charge features by summed intensity -- quick
            // sanity check that real proteoform charge fingerprints survived.
            if (multi > 0)
            {
                TestContext.Out.WriteLine("");
                TestContext.Out.WriteLine("Top 10 multi-charge features by summed intensity:");
                foreach (var f in features
                    .Where(f => f.ChargeCount >= 2)
                    .OrderByDescending(f => f.SummedIntensity)
                    .Take(10))
                {
                    string charges = string.Join(",", f.Charges.OrderBy(c => c));
                    TestContext.Out.WriteLine(
                        $"  mass={f.ConsensusMass,10:F4}  z={charges,-25}  rt={f.RTStart:F2}-{f.RTEnd:F2}  Sigma-intensity={f.SummedIntensity:E2}");
                }
            }
        }
    }
}
