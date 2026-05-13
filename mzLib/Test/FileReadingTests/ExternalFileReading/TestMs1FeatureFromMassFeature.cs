using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Deconvolution.Consensus;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.ExternalFileReading
{
    /// <summary>
    /// Round-trip coverage for <see cref="Ms1FeatureFile.FromMassFeatures"/>
    /// and the <see cref="Ms1FeatureFromMassFeatureExtensions.ToMs1Feature"/>
    /// extension. Builds fake <see cref="MassFeature"/> objects with known
    /// fields, writes them via the new factory, reads them back through the
    /// existing reader, and asserts every downstream-consumed field survives
    /// unchanged.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestMs1FeatureFromMassFeature
    {
        private static string directoryPath;

        [OneTimeSetUp]
        public void SetUp()
        {
            directoryPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\Ms1FeatureFromMassFeatureRoundTrip");
            Directory.CreateDirectory(directoryPath);
        }

        [OneTimeTearDown]
        public void TearDown()
        {
            Directory.Delete(directoryPath, true);
        }

        [Test]
        public static void ToMs1Feature_SingletonFeature_MapsAllFields()
        {
            // Single trace, single envelope. Apex == the one envelope.
            var feature = BuildMassFeature(id: 1, traces: new[]
            {
                BuildTrace(charge: 5, consensusMass: 1000.5,
                    envelopes: new[] { (rt: 10.5, intensity: 5.0e7) }),
            });

            var row = feature.ToMs1Feature(sequentialId: 0);

            Assert.That(row.SampleId, Is.EqualTo(0));
            Assert.That(row.Id, Is.EqualTo(0));
            Assert.That(row.Mass, Is.EqualTo(1000.5).Within(1e-9));
            Assert.That(row.Intensity, Is.EqualTo(5.0e7).Within(1e-3));
            Assert.That(row.RetentionTimeBegin, Is.EqualTo(10.5).Within(1e-9));
            Assert.That(row.RetentionTimeEnd, Is.EqualTo(10.5).Within(1e-9));
            Assert.That(row.RetentionTimeApex, Is.EqualTo(10.5).Within(1e-9));
            Assert.That(row.IntensityApex, Is.EqualTo(5.0e7).Within(1e-3));
            Assert.That(row.ChargeStateMin, Is.EqualTo(5));
            Assert.That(row.ChargeStateMax, Is.EqualTo(5));
            Assert.That(row.FractionIdMin, Is.EqualTo(0));
            Assert.That(row.FractionIdMax, Is.EqualTo(0));
        }

        [Test]
        public static void ToMs1Feature_MultiChargeFeature_ApexIsDominantTraceMaxEnvelope()
        {
            // Two traces; trace at charge=5 carries 2.4e8 total intensity (dominant),
            // trace at charge=6 carries 2.1e8. Apex must come from the charge=5 trace,
            // specifically its highest-intensity envelope at RT=15.2.
            var feature = BuildMassFeature(id: 2, traces: new[]
            {
                BuildTrace(charge: 5, consensusMass: 2000.6, envelopes: new[]
                {
                    (rt: 15.0, intensity: 8.0e7),
                    (rt: 15.2, intensity: 1.0e8), // apex of dominant trace
                    (rt: 15.4, intensity: 6.0e7),
                }),
                BuildTrace(charge: 6, consensusMass: 2000.8, envelopes: new[]
                {
                    (rt: 15.1, intensity: 1.2e8),
                    (rt: 15.3, intensity: 9.0e7),
                }),
            });

            var row = feature.ToMs1Feature(sequentialId: 1);

            // Cross-charge consensus is the intensity-weighted mean of per-trace
            // consensus masses: (2000.6 * 2.4e8 + 2000.8 * 2.1e8) / 4.5e8 ~= 2000.6933.
            Assert.That(row.Mass, Is.EqualTo((2000.6 * 2.4e8 + 2000.8 * 2.1e8) / 4.5e8).Within(1e-6));
            Assert.That(row.Intensity, Is.EqualTo(4.5e8).Within(1e-3));
            Assert.That(row.RetentionTimeBegin, Is.EqualTo(15.0).Within(1e-9));
            Assert.That(row.RetentionTimeEnd, Is.EqualTo(15.4).Within(1e-9));
            Assert.That(row.RetentionTimeApex, Is.EqualTo(15.2).Within(1e-9));
            Assert.That(row.IntensityApex, Is.EqualTo(1.0e8).Within(1e-3));
            Assert.That(row.ChargeStateMin, Is.EqualTo(5));
            Assert.That(row.ChargeStateMax, Is.EqualTo(6));
        }

        [Test]
        public static void ToMs1Feature_SampleIdAndFractionIdHonoured()
        {
            var feature = BuildMassFeature(id: 1, traces: new[]
            {
                BuildTrace(charge: 7, consensusMass: 500.25,
                    envelopes: new[] { (rt: 5.0, intensity: 1.0e6) }),
            });

            var row = feature.ToMs1Feature(sequentialId: 99, sampleId: 3, fractionId: 7);

            Assert.That(row.SampleId, Is.EqualTo(3));
            Assert.That(row.Id, Is.EqualTo(99));
            Assert.That(row.FractionIdMin, Is.EqualTo(7));
            Assert.That(row.FractionIdMax, Is.EqualTo(7));
        }

        [Test]
        public static void FromMassFeatures_WriteThenRead_AllFieldsSurvive()
        {
            var features = new[]
            {
                BuildMassFeature(id: 1, traces: new[]
                {
                    BuildTrace(charge: 5, consensusMass: 1000.5,
                        envelopes: new[] { (rt: 10.5, intensity: 5.0e7) }),
                }),
                BuildMassFeature(id: 2, traces: new[]
                {
                    BuildTrace(charge: 5, consensusMass: 2000.6, envelopes: new[]
                    {
                        (rt: 15.0, intensity: 8.0e7),
                        (rt: 15.2, intensity: 1.0e8),
                        (rt: 15.4, intensity: 6.0e7),
                    }),
                    BuildTrace(charge: 6, consensusMass: 2000.8, envelopes: new[]
                    {
                        (rt: 15.1, intensity: 1.2e8),
                        (rt: 15.3, intensity: 9.0e7),
                    }),
                }),
                BuildMassFeature(id: 3, traces: new[]
                {
                    BuildTrace(charge: 10, consensusMass: 5000.0, envelopes: new[]
                    {
                        (rt: 20.0, intensity: 2.0e7),
                        (rt: 20.5, intensity: 8.0e7),
                        (rt: 21.0, intensity: 3.0e7),
                    }),
                }),
            };

            string outputPath = Path.Combine(directoryPath, "consensus_ms1.feature");
            var written = Ms1FeatureFile.FromMassFeatures(features);
            written.WriteResults(outputPath);
            Assert.That(File.Exists(outputPath), Is.True);

            var roundTripped = FileReader.ReadFile<Ms1FeatureFile>(outputPath);
            Assert.That(roundTripped.Results.Count, Is.EqualTo(3));

            // Compare each record field-by-field against the original ToMs1Feature
            // output (the writer's canonical input). Tolerances are loose on
            // intensity and mass because CsvHelper round-trips doubles via
            // invariant-culture ToString/Parse -- some precision loss at the
            // last decimal is expected.
            for (int i = 0; i < features.Length; i++)
            {
                var original = features[i].ToMs1Feature(sequentialId: i);
                var read = roundTripped.Results[i];

                Assert.That(read.SampleId, Is.EqualTo(original.SampleId), $"row {i} SampleId");
                Assert.That(read.Id, Is.EqualTo(original.Id), $"row {i} Id");
                Assert.That(read.Mass, Is.EqualTo(original.Mass).Within(1e-9), $"row {i} Mass");
                Assert.That(read.Intensity, Is.EqualTo(original.Intensity).Within(1e-3), $"row {i} Intensity");
                Assert.That(read.RetentionTimeBegin, Is.EqualTo(original.RetentionTimeBegin).Within(1e-9), $"row {i} RtBegin");
                Assert.That(read.RetentionTimeEnd, Is.EqualTo(original.RetentionTimeEnd).Within(1e-9), $"row {i} RtEnd");
                Assert.That(read.RetentionTimeApex, Is.EqualTo(original.RetentionTimeApex).Within(1e-9), $"row {i} RtApex");
                Assert.That(read.IntensityApex, Is.EqualTo(original.IntensityApex).Within(1e-3), $"row {i} ApexIntensity");
                Assert.That(read.ChargeStateMin, Is.EqualTo(original.ChargeStateMin), $"row {i} ChargeMin");
                Assert.That(read.ChargeStateMax, Is.EqualTo(original.ChargeStateMax), $"row {i} ChargeMax");
                Assert.That(read.FractionIdMin, Is.EqualTo(original.FractionIdMin), $"row {i} FractionMin");
                Assert.That(read.FractionIdMax, Is.EqualTo(original.FractionIdMax), $"row {i} FractionMax");
            }
        }

        [Test]
        public static void FromMassFeatures_AssignsSequentialIdsFromZero()
        {
            var features = new[]
            {
                BuildMassFeature(id: 10, traces: new[] {
                    BuildTrace(charge: 5, consensusMass: 1000.0,
                        envelopes: new[] { (rt: 5.0, intensity: 1.0e6) }) }),
                BuildMassFeature(id: 20, traces: new[] {
                    BuildTrace(charge: 5, consensusMass: 2000.0,
                        envelopes: new[] { (rt: 6.0, intensity: 1.0e6) }) }),
                BuildMassFeature(id: 30, traces: new[] {
                    BuildTrace(charge: 5, consensusMass: 3000.0,
                        envelopes: new[] { (rt: 7.0, intensity: 1.0e6) }) }),
            };

            var file = Ms1FeatureFile.FromMassFeatures(features);

            // The factory ignores the producer's internal MassFeature.Id and
            // re-numbers from 0 to keep IDs dense and stable for downstream
            // consumers (a MassFeature pipeline may have non-contiguous IDs
            // after filtering, but the file should look clean).
            Assert.That(file.Results.Select(r => r.Id), Is.EqualTo(new[] { 0, 1, 2 }));
        }

        // ───────────────────────────────────────────────────────────────────
        // Fake-MassFeature builders
        // ───────────────────────────────────────────────────────────────────

        private static MassFeature BuildMassFeature(int id, IReadOnlyList<CorrectedTrace> traces)
        {
            var f = new MassFeature
            {
                Id = id,
                Traces = traces.ToList(),
            };
            f.Finalise();
            return f;
        }

        private static CorrectedTrace BuildTrace(
            int charge,
            double consensusMass,
            IReadOnlyList<(double rt, double intensity)> envelopes)
        {
            var trace = new CorrectedTrace
            {
                Id = charge, // arbitrary stable id for the test
                Charge = charge,
                ConsensusMass = consensusMass,
            };
            int scanIdx = 0;
            foreach (var (rt, intensity) in envelopes)
            {
                trace.Envelopes.Add(new CorrectedEnvelope
                {
                    ScanIndex = scanIdx,
                    ScanNumber = scanIdx + 1,
                    RT = rt,
                    OriginalMass = consensusMass,
                    CorrectedMass = consensusMass,
                    Charge = charge,
                    Intensity = intensity,
                    WasCorrected = false,
                });
                scanIdx++;
            }
            return trace;
        }
    }
}
