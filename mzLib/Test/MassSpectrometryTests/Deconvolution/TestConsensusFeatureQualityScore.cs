using Chemistry;
using MassSpectrometry;
using MassSpectrometry.Deconvolution.Consensus;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.MassSpectrometryTests.Deconvolution
{
    /// <summary>
    /// Tests for the per-feature deconvolution quality score carried from envelope, through the
    /// trace, onto <see cref="MassFeature"/>, and out to the <c>_ms1.feature</c> file.
    ///
    /// Groups:
    ///   A — the score reaches the feature at all, and NaN means "not scored"
    ///   B — the aggregate is intensity-weighted, and the maximum is reported alongside
    ///   C — the file round-trips, and the default written schema is unchanged
    /// </summary>
    [TestFixture]
    public class TestConsensusFeatureQualityScore
    {
        /// <summary>Two MS1 scans holding one species at one charge, so traces are predictable.</summary>
        private static (List<MsDataScan> Scans, List<IReadOnlyList<IsotopicEnvelope>> Envelopes)
            TwoScanFixture(double mass = 5000.0, int charge = 5,
                           double i1 = 1000, double i2 = 9000)
        {
            var scans = new List<MsDataScan>();
            var envs = new List<IReadOnlyList<IsotopicEnvelope>>();

            for (int k = 0; k < 2; k++)
            {
                double mz = mass.ToMz(charge);
                var spectrum = new MzSpectrum(
                    new[] { mz, mz + Constants.C13MinusC12 / charge },
                    new[] { 1e5, 5e4 }, false);
                scans.Add(new MsDataScan(spectrum, k + 1, 1, true, Polarity.Positive,
                    10.0 + k * 0.1, spectrum.Range, null, MZAnalyzerType.Orbitrap,
                    spectrum.SumOfAllY, null, null, null));
                envs.Add(new List<IsotopicEnvelope>
                {
                    new IsotopicEnvelope(mass, k == 0 ? i1 : i2, charge),
                });
            }
            return (scans, envs);
        }

        private static MassFeature FeatureFrom(
            List<MsDataScan> scans,
            List<IReadOnlyList<IsotopicEnvelope>> envs,
            Func<IsotopicEnvelope, MsDataScan, double> scorer)
        {
            var traces = MassTraceBuilder.BuildTraces(scans, envs, 0.05, 1, scorer);
            var corrected = traces.SelectMany(t => TraceCorrector.Correct(t)).ToList();
            var feature = new MassFeature { Id = 1 };
            feature.Traces.AddRange(corrected);
            feature.Finalise();
            return feature;
        }

        // ── A: the score arrives, and NaN means not scored ─────────────────────

        [Test]
        public void A1_NoScorer_LeavesScoresNaN()
        {
            var (scans, envs) = TwoScanFixture();
            var feature = FeatureFrom(scans, envs, null);

            Assert.That(double.IsNaN(feature.QualityScore), Is.True,
                "a feature traced without a scorer must report NaN, not 0");
            Assert.That(double.IsNaN(feature.MaxEnvelopeScore), Is.True);
            Assert.That(feature.Traces.SelectMany(t => t.Envelopes).All(e => double.IsNaN(e.Score)),
                Is.True, "every envelope should also be NaN");
        }

        [Test]
        public void A2_ScorerIsCalledOncePerEnvelope_AndValueReachesTheFeature()
        {
            var (scans, envs) = TwoScanFixture();
            int calls = 0;
            var feature = FeatureFrom(scans, envs, (e, s) => { calls++; return 0.75; });

            Assert.That(calls, Is.EqualTo(2), "one call per envelope observation");
            Assert.That(feature.QualityScore, Is.EqualTo(0.75).Within(1e-12));
            Assert.That(feature.MaxEnvelopeScore, Is.EqualTo(0.75).Within(1e-12));
        }

        [Test]
        public void A3_SpectrumAwareScorer_ProducesAScoreInRange()
        {
            var (scans, envs) = TwoScanFixture();
            var feature = FeatureFrom(scans, envs,
                DeconvolutionScorer.SpectrumAwareScorer(new Averagine()));

            Assert.That(double.IsNaN(feature.QualityScore), Is.False,
                "the real scorer should produce a value");
            Assert.That(feature.QualityScore, Is.InRange(0.0, 1.0));
        }

        // ── B: aggregation ────────────────────────────────────────────────────

        [Test]
        public void B1_AggregateIsIntensityWeighted_NotAPlainMean()
        {
            // Intensities 1000 and 9000; scores 0.0 and 1.0 respectively.
            // Plain mean would be 0.50; intensity-weighted is 9000/10000 = 0.90.
            var (scans, envs) = TwoScanFixture(i1: 1000, i2: 9000);
            var feature = FeatureFrom(scans, envs,
                (e, s) => e.TotalIntensity > 5000 ? 1.0 : 0.0);

            Assert.That(feature.QualityScore, Is.EqualTo(0.9).Within(1e-9),
                "a faint bad envelope must not outweigh an abundant good one");
            Assert.That(feature.QualityScore, Is.Not.EqualTo(0.5).Within(1e-3));
        }

        [Test]
        public void B2_MaximumIsReportedAlongsideTheMean()
        {
            var (scans, envs) = TwoScanFixture(i1: 1000, i2: 9000);
            var feature = FeatureFrom(scans, envs,
                (e, s) => e.TotalIntensity > 5000 ? 0.2 : 0.95);

            // One excellent faint envelope among abundant mediocre ones: the mean is dragged
            // down but the maximum preserves the fact that the species was seen cleanly once.
            Assert.That(feature.QualityScore, Is.LessThan(0.4));
            Assert.That(feature.MaxEnvelopeScore, Is.EqualTo(0.95).Within(1e-12));
        }

        [Test]
        public void B3_PartiallyScoredFeature_IgnoresTheUnscoredEnvelopes()
        {
            // A NaN must be excluded from the aggregate, not folded in as zero.
            var (scans, envs) = TwoScanFixture(i1: 1000, i2: 1000);
            var feature = FeatureFrom(scans, envs,
                (e, s) => s.OneBasedScanNumber == 1 ? 0.8 : double.NaN);

            Assert.That(feature.QualityScore, Is.EqualTo(0.8).Within(1e-12),
                "the unscored envelope should be skipped, not treated as 0");
            Assert.That(feature.MaxEnvelopeScore, Is.EqualTo(0.8).Within(1e-12));
        }

        // ── C: the file ───────────────────────────────────────────────────────

        [Test]
        public void C1_DefaultWrite_KeepsTheTopFdColumnSet()
        {
            var (scans, envs) = TwoScanFixture();
            var feature = FeatureFrom(scans, envs, (e, s) => 0.66);
            var file = Ms1FeatureFile.FromMassFeatures(new[] { feature });

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"qs_default_{Guid.NewGuid():N}_ms1.feature");
            try
            {
                file.WriteResults(path);
                string header = File.ReadLines(path).First();
                Assert.That(header, Does.Not.Contain("Quality_score"),
                    "the default schema must stay readable by tools expecting TopFD's columns");
                Assert.That(header, Does.Not.Contain("Max_envelope_score"));
                Assert.That(header, Does.Contain("Mass"));
            }
            finally { if (File.Exists(path)) File.Delete(path); }
        }

        [Test]
        public void C2_OptInWrite_EmitsTheScoreAndItRoundTrips()
        {
            var (scans, envs) = TwoScanFixture();
            var feature = FeatureFrom(scans, envs, (e, s) => 0.66);
            var file = Ms1FeatureFile.FromMassFeatures(new[] { feature });

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"qs_optin_{Guid.NewGuid():N}_ms1.feature");
            try
            {
                file.WriteResults(path, true);
                string header = File.ReadLines(path).First();
                Assert.That(header, Does.Contain("Quality_score"));

                var reread = new Ms1FeatureFile(path);
                var row = reread.Results.Single();
                Assert.That(row.QualityScore, Is.Not.Null);
                Assert.That(row.QualityScore.Value, Is.EqualTo(0.66).Within(1e-9));
                Assert.That(row.MaxEnvelopeScore.Value, Is.EqualTo(0.66).Within(1e-9));
            }
            finally { if (File.Exists(path)) File.Delete(path); }
        }

        [Test]
        public void C3_UnscoredFeature_WritesAnEmptyCellNotNaN()
        {
            // "NaN" in a numeric column breaks strict parsers; an empty cell is the honest
            // encoding of "this producer did not score its features".
            var (scans, envs) = TwoScanFixture();
            var feature = FeatureFrom(scans, envs, null);
            var file = Ms1FeatureFile.FromMassFeatures(new[] { feature });

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"qs_unscored_{Guid.NewGuid():N}_ms1.feature");
            try
            {
                file.WriteResults(path, true);
                string dataRow = File.ReadLines(path).Skip(1).First();
                Assert.That(dataRow, Does.Not.Contain("NaN"));

                var reread = new Ms1FeatureFile(path);
                Assert.That(reread.Results.Single().QualityScore, Is.Null);
            }
            finally { if (File.Exists(path)) File.Delete(path); }
        }

        [Test]
        public void C4_ExternalFileWithoutTheColumns_StillReads()
        {
            // A TopFD or FLASHDeconv file has no quality columns; reading one must not fail.
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"qs_external_{Guid.NewGuid():N}_ms1.feature");
            try
            {
                File.WriteAllLines(path, new[]
                {
                    "Sample_ID\tID\tMass\tIntensity\tTime_begin\tTime_end\tTime_apex\tApex_intensity\tMinimum_charge_state\tMaximum_charge_state\tMinimum_fraction_id\tMaximum_fraction_id",
                    "0\t1\t5000.0\t1000\t10.0\t10.2\t10.1\t900\t4\t6\t0\t0",
                });

                var reread = new Ms1FeatureFile(path);
                var row = reread.Results.Single();
                Assert.That(row.Mass, Is.EqualTo(5000.0).Within(1e-9));
                Assert.That(row.QualityScore, Is.Null, "absent column reads as null, not 0");
            }
            finally { if (File.Exists(path)) File.Delete(path); }
        }
    }
}
