using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using static Test.Deconvolution.DeconvolutionTestHelpers;

namespace Test.Deconvolution
{
    /// <summary>
    /// Tests that the deconvolution scorer handles negative charge states
    /// (ESI negative mode) correctly and produces equivalent results to
    /// the corresponding positive charge states.
    /// </summary>
    [TestFixture]
    public sealed class TestDeconvolutionScorerNegativeCharge
    {
        [Test]
        public void ComputeFeatures_NegativeCharge_CosineMatchesPositive()
        {
            // REGRESSION: failed before the fix
            var posEnv = BuildPerfectEnvelope(charge: TestCharge);
            var negEnv = BuildPerfectEnvelope(charge: -TestCharge);

            var posFeatures = DeconvolutionScorer.ComputeFeatures(posEnv, Model);
            var negFeatures = DeconvolutionScorer.ComputeFeatures(negEnv, Model);

            Assert.That(negFeatures.AveragineCosineSimilarity,
                Is.EqualTo(posFeatures.AveragineCosineSimilarity).Within(0.01),
                "Negative-charge cosine similarity must match positive-charge equivalent");
        }

        [Test]
        public void ComputeFeatures_NegativeCharge_PpmErrorMatchesPositive()
        {
            var posEnv = BuildPerfectEnvelope(charge: TestCharge);
            var negEnv = BuildPerfectEnvelope(charge: -TestCharge);

            var posFeatures = DeconvolutionScorer.ComputeFeatures(posEnv, Model);
            var negFeatures = DeconvolutionScorer.ComputeFeatures(negEnv, Model);

            Assert.That(negFeatures.AvgPpmError,
                Is.EqualTo(posFeatures.AvgPpmError).Within(0.1),
                "Negative-charge ppm error must match positive-charge equivalent");
        }

        [Test]
        public void ComputeFeatures_NegativeCharge_CompletenessMatchesPositive()
        {
            var posEnv = BuildPerfectEnvelope(charge: TestCharge);
            var negEnv = BuildPerfectEnvelope(charge: -TestCharge);

            var posFeatures = DeconvolutionScorer.ComputeFeatures(posEnv, Model);
            var negFeatures = DeconvolutionScorer.ComputeFeatures(negEnv, Model);

            Assert.That(negFeatures.PeakCompleteness,
                Is.EqualTo(posFeatures.PeakCompleteness).Within(0.01),
                "Negative-charge completeness must match positive-charge equivalent");
        }

        [Test]
        public void ComputeFeatures_NegativeCharge_RatioConsistencyMatchesPositive()
        {
            var posEnv = BuildPerfectEnvelope(charge: TestCharge);
            var negEnv = BuildPerfectEnvelope(charge: -TestCharge);

            var posFeatures = DeconvolutionScorer.ComputeFeatures(posEnv, Model);
            var negFeatures = DeconvolutionScorer.ComputeFeatures(negEnv, Model);

            Assert.That(negFeatures.IntensityRatioConsistency,
                Is.EqualTo(posFeatures.IntensityRatioConsistency).Within(0.01),
                "Negative-charge ratio consistency must match positive-charge equivalent");
        }

        [Test]
        public void ComputeScore_NegativeCharge_ScoreMatchesPositive()
        {
            var posEnv = BuildPerfectEnvelope(charge: TestCharge);
            var negEnv = BuildPerfectEnvelope(charge: -TestCharge);

            var posFeatures = DeconvolutionScorer.ComputeFeatures(posEnv, Model);
            var negFeatures = DeconvolutionScorer.ComputeFeatures(negEnv, Model);

            double posScore = DeconvolutionScorer.ComputeScore(posFeatures);
            double negScore = DeconvolutionScorer.ComputeScore(negFeatures);

            Assert.That(negScore, Is.EqualTo(posScore).Within(0.01),
                $"Negative-charge score ({negScore:F4}) must equal positive-charge score ({posScore:F4})");
        }

        [TestCase(-1)]
        [TestCase(-2)]
        [TestCase(-3)]
        [TestCase(-5)]
        [TestCase(-10)]
        public void ComputeFeatures_VariousNegativeCharges_AllFeaturesFinite(int negCharge)
        {
            var env = BuildPerfectEnvelope(charge: negCharge);
            var features = DeconvolutionScorer.ComputeFeatures(env, Model);

            Assert.That(double.IsFinite(features.AveragineCosineSimilarity), Is.True,
                $"Cosine must be finite for charge={negCharge}");
            Assert.That(double.IsFinite(features.AvgPpmError), Is.True,
                $"PpmError must be finite for charge={negCharge}");
            Assert.That(double.IsFinite(features.PeakCompleteness), Is.True,
                $"Completeness must be finite for charge={negCharge}");
            Assert.That(double.IsFinite(features.IntensityRatioConsistency), Is.True,
                $"RatioConsistency must be finite for charge={negCharge}");
        }

        [TestCase(-1)]
        [TestCase(-3)]
        [TestCase(-5)]
        [TestCase(-10)]
        public void ComputeFeatures_NegativeCharge_CosineAboveThreshold(int negCharge)
        {
            var env = BuildPerfectEnvelope(charge: negCharge);
            var features = DeconvolutionScorer.ComputeFeatures(env, Model);

            Assert.That(features.AveragineCosineSimilarity, Is.GreaterThanOrEqualTo(0.99),
                $"Perfect negative-charge envelope (z={negCharge}) should have cosine >= 0.99. Got {features.AveragineCosineSimilarity:F4}");
        }

        [Test]
        public void ScoreEnvelopes_MixedPolarity_ReturnsCorrectCount()
        {
            var envelopes = new[]
            {
                BuildPerfectEnvelope(TestMass, TestCharge),
                BuildPerfectEnvelope(TestMass, -TestCharge),
                BuildPerfectEnvelope(TestMass + 1000.0, -3),
            };

            var pairs = DeconvolutionScorer.ScoreEnvelopes(envelopes, Model).ToList();

            Assert.That(pairs.Count, Is.EqualTo(3),
                "ScoreEnvelopes should return one pair per envelope regardless of charge sign");

            for (int i = 0; i < envelopes.Length; i++)
            {
                Assert.That(pairs[i].Envelope, Is.SameAs(envelopes[i]),
                    $"Pair {i} must reference the original envelope object");
            }
        }

        [Test]
        public void ScoreEnvelopes_NegativeChargeEnvelope_ScoreAboveHalf()
        {
            var env = BuildPerfectEnvelope(charge: -TestCharge);

            var pairs = DeconvolutionScorer.ScoreEnvelopes(new[] { env }, Model).ToList();

            Assert.That(pairs[0].Score, Is.GreaterThan(0.5),
                $"Perfect negative-charge envelope should score above 0.5. Got {pairs[0].Score:F4}");
        }

        [Test]
        public void BuildPerfectEnvelope_NegativeCharge_PeaksAscendInMz()
        {
            // REGRESSION: failed before the fix — isotope step was negative
            var env = BuildPerfectEnvelope(charge: -TestCharge);

            for (int i = 1; i < env.Peaks.Count; i++)
            {
                Assert.That(env.Peaks[i].mz, Is.GreaterThan(env.Peaks[i - 1].mz),
                    $"Peak {i} m/z ({env.Peaks[i].mz:F6}) must be greater than peak {i - 1} m/z ({env.Peaks[i - 1].mz:F6}) for negative charge envelope");
            }
        }
    }
}
