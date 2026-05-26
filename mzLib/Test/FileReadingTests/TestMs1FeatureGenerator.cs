using System.IO;
using System.Linq;
using MassSpectrometry.Deconvolution.Consensus;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Tests the write side of the whole-file driver (<see cref="Ms1FeatureGenerator"/>):
    /// retention time is emitted in SECONDS (to match real OpenMS FLASHDeconv output) and the
    /// header is the TopFD-superset <c>_ms1.feature</c> schema (the 11 canonical FLASHDeconv
    /// columns plus Apex_intensity), round-tripping cleanly through <see cref="Ms1FeatureFile"/>.
    /// The deconvolve+trace half of the driver is covered end-to-end by the [Explicit]
    /// JurkatMetaFlashDeconGenerator (needs local data).
    /// </summary>
    [TestFixture]
    public class TestMs1FeatureGenerator
    {
        // CsvHelper writes every Ms1Feature property in declaration order (first [Name] alias):
        // the 11 canonical FLASHDeconv columns plus Apex_intensity (the TopFD superset).
        private const string ExpectedHeader =
            "Sample_ID\tID\tMass\tIntensity\tTime_begin\tTime_end\tTime_apex\tApex_intensity\t" +
            "Minimum_charge_state\tMaximum_charge_state\tMinimum_fraction_id\tMaximum_fraction_id";

        // The columns the REAL OpenMS FLASHDeconv _ms1.feature carries (no Apex_intensity).
        private static readonly string[] RealFlashDeconvColumns =
        {
            "Sample_ID", "ID", "Mass", "Intensity", "Time_begin", "Time_end", "Time_apex",
            "Minimum_charge_state", "Maximum_charge_state", "Minimum_fraction_id", "Maximum_fraction_id"
        };

        [Test]
        public void WriteFeatures_EmitsSecondsRt_TopFdSupersetSchema_RoundTrips()
        {
            // One feature: neutral mass 5000 at charge 5, two scans at RT 10.0 and 10.2 MINUTES.
            var feature = SyntheticFeature(mass: 5000.0, charge: 5, rtMinutes: new[] { 10.0, 10.2 }, intensity: 1e7);
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "metaflashdecon_roundtrip_ms1.feature");
            try
            {
                Ms1FeatureGenerator.WriteFeatures(new[] { feature }, path, secondsPerRtUnit: 60.0);

                // Header is the exact TopFD-superset schema and contains every real FLASHDeconv column.
                string header = File.ReadLines(path).First();
                Assert.That(header, Is.EqualTo(ExpectedHeader));
                var cols = header.Split('\t');
                foreach (var col in RealFlashDeconvColumns)
                    Assert.That(cols, Contains.Item(col), $"missing FLASHDeconv column {col}");

                // Round-trip the values; RT must be in seconds (minutes x 60).
                var roundTrip = FileReader.ReadFile<Ms1FeatureFile>(path);
                Assert.That(roundTrip.Results, Has.Count.EqualTo(1));
                var row = roundTrip.Results[0];
                Assert.That(row.Mass, Is.EqualTo(5000.0).Within(1e-6));
                Assert.That(row.ChargeStateMin, Is.EqualTo(5));
                Assert.That(row.ChargeStateMax, Is.EqualTo(5));
                Assert.That(row.RetentionTimeBegin, Is.EqualTo(600.0).Within(1e-6), "10.0 min -> 600 s");
                Assert.That(row.RetentionTimeEnd, Is.EqualTo(612.0).Within(1e-6), "10.2 min -> 612 s");
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [Test]
        public void WriteFeatures_ScaleOne_LeavesRtInSourceUnits()
        {
            var feature = SyntheticFeature(5000.0, 5, new[] { 10.0, 10.2 }, 1e7);
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "metaflashdecon_minutes_ms1.feature");
            try
            {
                Ms1FeatureGenerator.WriteFeatures(new[] { feature }, path, secondsPerRtUnit: 1.0);
                var row = FileReader.ReadFile<Ms1FeatureFile>(path).Results.Single();
                Assert.That(row.RetentionTimeBegin, Is.EqualTo(10.0).Within(1e-6));
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        private static MassFeature SyntheticFeature(double mass, int charge, double[] rtMinutes, double intensity)
        {
            var envs = rtMinutes.Select((rt, i) => new CorrectedEnvelope
            {
                ScanIndex = i,
                ScanNumber = i + 1,
                RT = rt,
                OriginalMass = mass,
                CorrectedMass = mass,
                Charge = charge,
                Intensity = intensity,
                WasCorrected = false,
            }).ToList();

            var feature = new MassFeature { Id = 1 };
            feature.Traces.Add(new CorrectedTrace { Id = charge, Charge = charge, ConsensusMass = mass, Envelopes = envs });
            feature.Finalise();
            return feature;
        }
    }
}
