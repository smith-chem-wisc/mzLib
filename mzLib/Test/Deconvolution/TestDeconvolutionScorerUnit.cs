using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using static Test.Deconvolution.DeconvolutionTestHelpers;

namespace Test.Deconvolution
{
    /// <summary>
    /// Unit tests for <see cref="DeconvolutionScorer"/> and <see cref="EnvelopeScoreFeatures"/>.
    /// These tests operate on synthetic envelopes built directly from the Averagine model
    /// — no deconvolution algorithm is called. Tests are grouped by feature (A), score (B),
    /// and the ScoreEnvelopes API (C).
    /// </summary>
    [TestFixture]
    public sealed class TestDeconvolutionScorerUnit
    {
        // ══════════════════════════════════════════════════════════════════════
        // Group A — ComputeFeatures: feature correctness on synthetic data
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void A1_ComputeFeatures_PerfectSyntheticEnvelope_CosineNearOne()
        {
            var env      = BuildPerfectEnvelope();
            var features = DeconvolutionScorer.ComputeFeatures(env, Model);

            Assert.That(features.AveragineCosineSimilarity, Is.GreaterThanOrEqualTo(0.99),
                $"Perfect synthetic envelope should have cosine ≥ 0.99. Got {features.AveragineCosineSimilarity:F4}");

            Assert.That(features.IntensityRatioConsistency, Is.GreaterThanOrEqualTo(0.90),
                $"Perfect synthetic envelope should have IntensityRatioConsistency ≥ 0.90. Got {features.IntensityRatioConsistency:F4}");
        }

        [Test]
        public void A2_ComputeFeatures_PerfectSyntheticEnvelope_PpmErrorNearZero()
        {
            var env      = BuildPerfectEnvelope();
            var features = DeconvolutionScorer.ComputeFeatures(env, Model);

            Assert.That(features.AvgPpmError, Is.LessThan(0.01),
                $"Peaks at exact theoretical positions should have AvgPpmError < 0.01 ppm. Got {features.AvgPpmError:F4}");
        }

        [Test]
        public void A3_ComputeFeatures_PerfectSyntheticEnvelope_CompletenessIsOne()
        {
            var env      = BuildPerfectEnvelope();
            var features = DeconvolutionScorer.ComputeFeatures(env, Model);

            Assert.That(features.PeakCompleteness, Is.EqualTo(1.0).Within(0.01),
                $"All expected peaks present: completeness should be 1.0. Got {features.PeakCompleteness:F4}");
        }

        [Test]
        public void A4_ComputeFeatures_MissingHalfPeaks_CompletenessIsHalfApprox()
        {
            // Build an envelope retaining only even-index isotope peaks (0, 2, 4, …)
            var full = BuildPerfectEnvelope();

            var halfPeaks = full.Peaks
                .Select((p, i) => (p, i))
                .Where(x => x.i % 2 == 0)
                .Select(x => x.p)
                .ToList();

            var env = new IsotopicEnvelope(
                id:               0,
                peaks:            halfPeaks,
                monoisotopicmass: full.MonoisotopicMass,
                chargestate:      full.Charge,
                intensity:        halfPeaks.Sum(p => p.intensity),
                score:            0.0);

            var features = DeconvolutionScorer.ComputeFeatures(env, Model);

            Assert.That(features.PeakCompleteness, Is.InRange(0.3, 0.7),
                $"Half the peaks present: completeness should be near 0.5. Got {features.PeakCompleteness:F4}");
        }

        [Test]
        public void A5_ComputeFeatures_ShiftedPeaks_CosineDrops()
        {
            // Shift every peak m/z by 50 ppm — far outside the 10 ppm matching window.
            // The observed vector will be all zeros and cosine should be 0.
            var full = BuildPerfectEnvelope();

            var shiftedPeaks = full.Peaks
                .Select(p => (p.mz * (1.0 + 50e-6), p.intensity))
                .ToList();

            var env = new IsotopicEnvelope(
                id:               0,
                peaks:            shiftedPeaks,
                monoisotopicmass: full.MonoisotopicMass,
                chargestate:      full.Charge,
                intensity:        full.TotalIntensity,
                score:            0.0);

            var features = DeconvolutionScorer.ComputeFeatures(env, Model);

            Assert.That(features.AveragineCosineSimilarity, Is.LessThan(0.5),
                $"50 ppm shifted peaks should not match Averagine positions. Cosine: {features.AveragineCosineSimilarity:F4}");
        }

        [Test]
        public void A6_ComputeFeatures_SingleZeroIntensityPeak_DoesNotThrow()
        {
            // A degenerate envelope with a single zero-intensity peak
            var peaks = new List<(double mz, double intensity)>
                { (TestMass.ToMz(TestCharge), 0.0) };

            var env = new IsotopicEnvelope(
                id:               0,
                peaks:            peaks,
                monoisotopicmass: TestMass,
                chargestate:      TestCharge,
                intensity:        0.0,
                score:            0.0);

            EnvelopeScoreFeatures features = default;
            Assert.DoesNotThrow(() => features = DeconvolutionScorer.ComputeFeatures(env, Model),
                "ComputeFeatures must not throw for a degenerate single-peak envelope");

            Assert.That(double.IsFinite(features.AveragineCosineSimilarity), Is.True);
            Assert.That(double.IsFinite(features.AvgPpmError),               Is.True);
            Assert.That(double.IsFinite(features.PeakCompleteness),          Is.True);
            Assert.That(double.IsFinite(features.IntensityRatioConsistency), Is.True);
        }

        // ══════════════════════════════════════════════════════════════════════
        // Group B — ComputeScore: score ordering and bounds
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void B1_ComputeScore_HighQualityFeatures_ScoreAbovePointFive()
        {
            var features = new EnvelopeScoreFeatures(
                averagineCosineSimilarity: 0.95,
                avgPpmError:              1.5,
                peakCompleteness:         0.95,
                intensityRatioConsistency: 0.95);

            double score = DeconvolutionScorer.ComputeScore(features);

            Assert.That(score, Is.GreaterThan(0.5),
                $"High-quality features should produce score > 0.5. Got {score:F4}");
        }

        [Test]
        public void B2_ComputeScore_LowQualityFeatures_ScoreBelowPointFive()
        {
            var features = new EnvelopeScoreFeatures(
                averagineCosineSimilarity: 0.2,
                avgPpmError:              25.0,
                peakCompleteness:         0.2,
                intensityRatioConsistency: 0.1);

            double score = DeconvolutionScorer.ComputeScore(features);

            Assert.That(score, Is.LessThan(0.5),
                $"Low-quality features should produce score < 0.5. Got {score:F4}");
        }

        [Test]
        public void B3_ComputeScore_ScoreIsInZeroOneRange(
            [Values(0.0, 0.5, 1.0)] double cosine,
            [Values(0.0, 10.0, 100.0)] double ppm,
            [Values(0.0, 0.5, 1.0)] double completeness,
            [Values(0.0, 0.5, 1.0)] double ratioConsistency)
        {
            var features = new EnvelopeScoreFeatures(cosine, ppm, completeness, ratioConsistency);
            double score = DeconvolutionScorer.ComputeScore(features);

            Assert.That(score, Is.InRange(0.0, 1.0),
                $"Score must be in [0,1] for cosine={cosine}, ppm={ppm}, completeness={completeness}, ratioConsistency={ratioConsistency}. Got {score:F6}");
        }

        [Test]
        public void B4_ComputeScore_BetterFeaturesGiveHigherScore()
        {
            var highQ = new EnvelopeScoreFeatures(0.95, 1.0, 0.95, 0.95);
            var lowQ  = new EnvelopeScoreFeatures(0.30, 20.0, 0.30, 0.15);

            double scoreHigh = DeconvolutionScorer.ComputeScore(highQ);
            double scoreLow  = DeconvolutionScorer.ComputeScore(lowQ);

            Assert.That(scoreHigh, Is.GreaterThan(scoreLow),
                $"High-quality features ({scoreHigh:F4}) should score above low-quality ({scoreLow:F4})");
        }

        // ══════════════════════════════════════════════════════════════════════
        // Group C — ScoreEnvelopes: API correctness
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void C1_ScoreEnvelopes_ReturnsPairPerEnvelope()
        {
            var envelopes = new[]
            {
                BuildPerfectEnvelope(TestMass,           TestCharge),
                BuildPerfectEnvelope(TestMass + 1000.0,  TestCharge + 1),
                BuildPerfectEnvelope(TestMass + 2000.0,  TestCharge + 2),
            };

            var pairs = DeconvolutionScorer.ScoreEnvelopes(envelopes, Model).ToList();

            Assert.That(pairs.Count, Is.EqualTo(3),
                "ScoreEnvelopes should return one pair per input envelope");

            for (int i = 0; i < envelopes.Length; i++)
                Assert.That(pairs[i].Envelope, Is.SameAs(envelopes[i]),
                    $"Pair {i} Envelope reference should be the original object");
        }

        [Test]
        public void C2_ScoreEnvelopes_PreservesOriginalEnvelopeScore()
        {
            var env = BuildPerfectEnvelope();
            double originalScore = env.Score;

            // Score via the generic scorer
            _ = DeconvolutionScorer.ScoreEnvelopes(new[] { env }, Model).ToList();

            Assert.That(env.Score, Is.EqualTo(originalScore),
                "ScoreEnvelopes must not mutate the original envelope's Score field");
        }

        [Test]
        public void C3_ScoreEnvelopes_EmptyInput_ReturnsEmpty()
        {
            var result = DeconvolutionScorer.ScoreEnvelopes(
                Enumerable.Empty<IsotopicEnvelope>(), Model).ToList();

            Assert.That(result, Is.Empty,
                "ScoreEnvelopes on an empty sequence must return an empty sequence");
        }

        // ══════════════════════════════════════════════════════════════════════
        // Group D — DeconvoluteWithGenericScoring: integration tests
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void D1_DeconvoluteWithGenericScoring_MzSpectrum_SetsGenericScoreOnEveryEnvelope()
        {
            // Arrange: build a spectrum from a perfect synthetic envelope at charge 2.
            // Charge 2 needs 2/5 = 0 adjacent charge states observed by Classic deconvolution,
            // so it works with an isolated synthetic isotope pattern.
            var env = BuildPerfectEnvelope(monoMass: TestMass, charge: 2);
            double[] mzArray = env.Peaks.Select(p => p.mz).ToArray();
            double[] intensityArray = env.Peaks.Select(p => p.intensity).ToArray();
            var spectrum = new MzSpectrum(mzArray, intensityArray, false);

            var deconParams = new ClassicDeconvolutionParameters(
                minCharge: 1, maxCharge: 10, deconPpm: 20, intensityRatio: 3);

            // Act
            var results = Deconvoluter.DeconvoluteWithGenericScoring(spectrum, deconParams).ToList();

            // Assert
            Assert.That(results, Is.Not.Empty, "Should deconvolute at least one envelope");
            foreach (var envelope in results)
            {
                Assert.That(envelope.HasGenericScore, Is.True,
                    "Every yielded envelope must have GenericScore set");
                Assert.That(envelope.GenericScore!.Value, Is.GreaterThan(0.5),
                    $"Perfect synthetic Averagine envelope must score > 0.5. Got {envelope.GenericScore}");
            }
        }

        [Test]
        public void D2_DeconvoluteWithGenericScoring_MzSpectrum_DoesNotMutateAlgorithmScore()
        {
            // Arrange: same spectrum/params as D1
            var env = BuildPerfectEnvelope(monoMass: TestMass, charge: 2);
            double[] mzArray = env.Peaks.Select(p => p.mz).ToArray();
            double[] intensityArray = env.Peaks.Select(p => p.intensity).ToArray();

            var deconParams = new ClassicDeconvolutionParameters(
                minCharge: 1, maxCharge: 10, deconPpm: 20, intensityRatio: 3);

            // Get baseline scores from plain Deconvolute
            var spectrum1 = new MzSpectrum(mzArray, intensityArray, false);
            var baseline = Deconvoluter.Deconvolute(spectrum1, deconParams).ToList();

            // Get scored results from a fresh spectrum
            var spectrum2 = new MzSpectrum(mzArray, intensityArray, false);
            var scored = Deconvoluter.DeconvoluteWithGenericScoring(spectrum2, deconParams).ToList();

            // Assert: same number of envelopes, and Score values match
            Assert.That(scored.Count, Is.EqualTo(baseline.Count),
                "DeconvoluteWithGenericScoring should yield the same number of envelopes as Deconvolute");
            for (int i = 0; i < baseline.Count; i++)
            {
                Assert.That(scored[i].Score, Is.EqualTo(baseline[i].Score).Within(1e-10),
                    $"Envelope {i}: algorithm-specific Score must be unchanged by generic scoring");
            }
        }

        [Test]
        public void D3_DeconvoluteWithGenericScoring_MsDataScan_SetsGenericScore()
        {
            // Arrange: same envelope as D1, wrapped in an MsDataScan
            var env = BuildPerfectEnvelope(monoMass: TestMass, charge: 2);
            double[] mzArray = env.Peaks.Select(p => p.mz).ToArray();
            double[] intensityArray = env.Peaks.Select(p => p.intensity).ToArray();
            var spectrum = new MzSpectrum(mzArray, intensityArray, false);

            var scan = new MsDataScan(spectrum, 1, 1, true, Polarity.Positive, 1.0,
                spectrum.Range, "test scan", MZAnalyzerType.Unknown, spectrum.SumOfAllY,
                null, null, null);

            var deconParams = new ClassicDeconvolutionParameters(
                minCharge: 1, maxCharge: 10, deconPpm: 20, intensityRatio: 3);

            // Act
            var results = Deconvoluter.DeconvoluteWithGenericScoring(scan, deconParams).ToList();

            // Assert
            Assert.That(results, Is.Not.Empty, "MsDataScan overload should yield envelopes");
            foreach (var envelope in results)
            {
                Assert.That(envelope.HasGenericScore, Is.True,
                    "MsDataScan overload must also set GenericScore on every envelope");
                Assert.That(envelope.GenericScore!.Value, Is.GreaterThan(0.5),
                    $"Perfect synthetic Averagine envelope must score > 0.5. Got {envelope.GenericScore}");
            }
        }
    }
}
