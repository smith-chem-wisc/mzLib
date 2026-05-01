using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using static Test.MassSpectrometryTests.Deconvolution.DeconvolutionTestHelpers;

namespace Test.MassSpectrometryTests.Deconvolution
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

            // BuildPerfectEnvelope scales every peak by the same baseIntens, so the
            // per-peak ratio CV is exactly 0 and IntensityRatioConsistency = 1/(1+0) = 1.0.
            // The construction pins this mathematically; the assertion should pin it
            // tightly so a regression introducing even small per-peak ratio drift fails.
            Assert.That(features.IntensityRatioConsistency, Is.EqualTo(1.0).Within(1e-6),
                $"Perfect synthetic envelope must produce IntensityRatioConsistency == 1.0. Got {features.IntensityRatioConsistency:F6}");
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

            // The 10 ppm matching window means a 50 ppm shift produces an all-zero
            // observed vector; cosine must be exactly 0. A looser bound (e.g.
            // < 0.5) would silently accept any partial leakage of out-of-window
            // intensity into observed[], so pin the contract precisely.
            Assert.That(features.AveragineCosineSimilarity, Is.EqualTo(0.0).Within(1e-9),
                $"50 ppm shifted peaks must produce cosine == 0. Got: {features.AveragineCosineSimilarity:F6}");
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
        // Group D — UseGenericScore via Deconvolute: integration tests
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void D1_Deconvolute_UseGenericScore_MzSpectrum_SetsGenericScoreOnEveryEnvelope()
        {
            // Arrange: build a spectrum from a perfect synthetic envelope at charge 2.
            // Charge 2 needs 2/5 = 0 adjacent charge states observed by Classic deconvolution,
            // so it works with an isolated synthetic isotope pattern.
            var env = BuildPerfectEnvelope(monoMass: TestMass, charge: 2);
            double[] mzArray = env.Peaks.Select(p => p.mz).ToArray();
            double[] intensityArray = env.Peaks.Select(p => p.intensity).ToArray();
            var spectrum = new MzSpectrum(mzArray, intensityArray, false);

            var deconParams = new ClassicDeconvolutionParameters(
                minCharge: 1, maxCharge: 10, deconPpm: 20, intensityRatio: 3)
            { UseGenericScore = true };

            // Act
            var results = Deconvoluter.Deconvolute(spectrum, deconParams).ToList();

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
        public void D2_Deconvolute_UseGenericScore_MzSpectrum_DoesNotMutateAlgorithmScore()
        {
            // Arrange: same spectrum as D1, but now with two parameter objects --
            // one with UseGenericScore off (baseline) and one with it on (scored).
            var env = BuildPerfectEnvelope(monoMass: TestMass, charge: 2);
            double[] mzArray = env.Peaks.Select(p => p.mz).ToArray();
            double[] intensityArray = env.Peaks.Select(p => p.intensity).ToArray();

            var baselineParams = new ClassicDeconvolutionParameters(
                minCharge: 1, maxCharge: 10, deconPpm: 20, intensityRatio: 3);
            var scoringParams = new ClassicDeconvolutionParameters(
                minCharge: 1, maxCharge: 10, deconPpm: 20, intensityRatio: 3)
            { UseGenericScore = true };

            // Get baseline scores from plain Deconvolute (UseGenericScore = false)
            var spectrum1 = new MzSpectrum(mzArray, intensityArray, false);
            var baseline = Deconvoluter.Deconvolute(spectrum1, baselineParams).ToList();

            // Get scored results from a fresh spectrum with UseGenericScore = true
            var spectrum2 = new MzSpectrum(mzArray, intensityArray, false);
            var scored = Deconvoluter.Deconvolute(spectrum2, scoringParams).ToList();

            // Assert: same number of envelopes, and algorithm-specific Score values match
            Assert.That(scored.Count, Is.EqualTo(baseline.Count),
                "UseGenericScore must not change the number of envelopes yielded");
            for (int i = 0; i < baseline.Count; i++)
            {
                Assert.That(scored[i].Score, Is.EqualTo(baseline[i].Score).Within(1e-10),
                    $"Envelope {i}: algorithm-specific Score must be unchanged by generic scoring");
            }
        }

        [Test]
        public void D3_Deconvolute_UseGenericScore_MsDataScan_SetsGenericScore()
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
                minCharge: 1, maxCharge: 10, deconPpm: 20, intensityRatio: 3)
            { UseGenericScore = true };

            // Act
            var results = Deconvoluter.Deconvolute(scan, deconParams).ToList();

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
