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
    /// stitching) against a local raw/mzML and write the result as a FLASHDeconv-
    /// style <c>_ms1.feature</c> file via <see cref="Ms1FeatureFile.FromMassFeatures"/>.
    ///
    /// Lives as an <c>[Explicit]</c> test because each case consumes a hundreds-of-MB
    /// local raw/mzML at a hard-coded path that's not part of the repo or CI.
    /// Run with:
    /// <code>
    /// dotnet test --filter FullyQualifiedName~TestConsensusToMs1FeatureEndToEnd
    /// </code>
    ///
    /// Pipeline parameters match the values established in the consensus-tracing
    /// research arc: Pass B grouping at 1.5 Da (intentionally loose to cohort
    /// off-by-one twins for correction), 1-scan gap tolerance, 10 ppm cross-charge
    /// mass agreement. The Classic charge range is per-case because top-down and
    /// bottom-up have very different effective ranges.
    /// </summary>
    [TestFixture]
    [Explicit("End-to-end consensus pipeline driver; needs local raw/mzML at E:\\TestData\\MetaMorpheus\\.")]
    [ExcludeFromCodeCoverage]
    public class TestConsensusToMs1FeatureEndToEnd
    {
        private const double PassBToleranceDa = 1.5;
        private const int    MaxGap            = 1;
        private const double CrossChargeMassPpm = 10.0;

        /// <summary>
        /// One driver invocation against a single mzML/raw file. Parameters are
        /// captured per-case because top-down vs bottom-up wants very different
        /// Classic charge ranges (1..60 vs 1..12 here).
        /// </summary>
        public sealed record DriverSpec(string Label, string SourcePath, string OutputPath, int MaxCharge)
        {
            public override string ToString() => Label;
        }

        public static IEnumerable<DriverSpec> Specs()
        {
            // Bottom-up HEK293 on Velos. Tryptic peptides rarely go above +6,
            // so 12 is a comfortable ceiling matching MetaMorpheus's standard
            // DeconvolutionMaxAssumedChargeState = 12 default for BU runs.
            yield return new DriverSpec(
                "BU_293_Velos_3",
                @"E:\TestData\MetaMorpheus\BottomUp\20100609_Velos1_TaGe_SA_293_3.raw",
                @"E:\TestData\MetaMorpheus\BottomUp\20100609_Velos1_TaGe_SA_293_3_ms1.feature",
                MaxCharge: 12);
            yield return new DriverSpec(
                "BU_293_Velos_4",
                @"E:\TestData\MetaMorpheus\BottomUp\20100609_Velos1_TaGe_SA_293_4.raw",
                @"E:\TestData\MetaMorpheus\BottomUp\20100609_Velos1_TaGe_SA_293_4_ms1.feature",
                MaxCharge: 12);
        }

        [Test]
        [TestCaseSource(nameof(Specs))]
        public void RunFullPipelineAndEmitMs1FeatureFile(DriverSpec spec)
        {
            Assume.That(File.Exists(spec.SourcePath),
                $"Source not found at {spec.SourcePath}; install local test data or override the spec.");

            var sw = Stopwatch.StartNew();

            // 1. Load the file and extract MS1 scans in order.
            var dataFile = MsDataFileReader.GetDataFile(spec.SourcePath).LoadAllStaticData();
            var ms1Scans = dataFile.GetAllScansList()
                .Where(s => s.MsnOrder == 1)
                .OrderBy(s => s.OneBasedScanNumber)
                .ToList();
            Assume.That(ms1Scans, Is.Not.Empty);
            TestContext.Out.WriteLine($"[{spec.Label}] MS1 scans loaded: {ms1Scans.Count}  (elapsed {sw.Elapsed})");

            // 2. Per-scan Classic decon. Charge range chosen to match the
            //    MetaMorpheus task config used for the corresponding benchmark.
            var deconParams = new ClassicDeconvolutionParameters(
                minCharge: 1,
                maxCharge: spec.MaxCharge,
                deconPpm: 4.0,
                intensityRatio: 3.0);
            var perScanEnvelopes = ms1Scans
                .Select(s => Deconvoluter.Deconvolute(s.MassSpectrum, deconParams).ToList())
                .ToList();
            int totalEnvelopes = perScanEnvelopes.Sum(s => s.Count);
            TestContext.Out.WriteLine($"[{spec.Label}] Classic envelopes total: {totalEnvelopes}  (elapsed {sw.Elapsed})");

            // 3. Pass B trace grouping (loose: 1.5 Da) so off-by-one twins
            //    cohort with their parents and the corrector can rescue them.
            var traces = MassTraceBuilder.BuildTraces(
                ms1Scans,
                perScanEnvelopes.Select(scan => (IReadOnlyList<IsotopicEnvelope>)scan).ToList(),
                PassBToleranceDa,
                MaxGap);
            TestContext.Out.WriteLine($"[{spec.Label}] Traces built: {traces.Count}  (elapsed {sw.Elapsed})");

            // 4. Per-trace weighted-median off-by-one correction. Uniform
            //    weight matches NOTES.md's Phase 5 finding that weighting
            //    scheme barely changes consensus mass.
            var corrected = traces.Select(t => TraceCorrector.Correct(t)).ToList();
            int corrections = corrected.Sum(c => c.CorrectionCount);
            TestContext.Out.WriteLine($"[{spec.Label}] Corrections applied: {corrections}  (elapsed {sw.Elapsed})");

            // 5. Cross-charge clustering.
            var features = MassFeatureBuilder.BuildFeatures(corrected, CrossChargeMassPpm);
            int multi = features.Count(f => f.ChargeCount >= 2);
            TestContext.Out.WriteLine($"[{spec.Label}] MassFeatures: {features.Count} (multi-charge: {multi})  (elapsed {sw.Elapsed})");

            // 6. Write FLASHDeconv-style _ms1.feature via the new factory.
            Ms1FeatureFile.FromMassFeatures(features).WriteResults(spec.OutputPath);

            sw.Stop();
            var bytes = new FileInfo(spec.OutputPath).Length;
            TestContext.Out.WriteLine($"[{spec.Label}] Wrote {bytes / 1024.0:F1} KB to {spec.OutputPath}  (total elapsed {sw.Elapsed})");

            // Sanity: the file should round-trip cleanly.
            var roundTrip = FileReader.ReadFile<Ms1FeatureFile>(spec.OutputPath);
            Assert.That(roundTrip.Results.Count, Is.EqualTo(features.Count),
                "Round-trip row count mismatch -- writer or reader regression.");
            TestContext.Out.WriteLine($"[{spec.Label}] Round-trip read: {roundTrip.Results.Count} rows.");

            // Top 10 multi-charge features by summed intensity -- quick
            // sanity check that real proteoform/peptide charge fingerprints survived.
            if (multi > 0)
            {
                TestContext.Out.WriteLine("");
                TestContext.Out.WriteLine($"[{spec.Label}] Top 10 multi-charge features by summed intensity:");
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
