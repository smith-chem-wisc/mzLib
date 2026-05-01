using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.MassSpectrometryTests.Deconvolution
{
    /// <summary>
    /// Tests that validate DeconvolutionScorer against independently sourced
    /// isotope envelopes, ensuring the scorer produces physically correct results
    /// rather than merely agreeing with its own Averagine model.
    ///
    /// The DeconvolutionScorer's logistic weights were derived by fitting several
    /// hundred thousand top-down envelopes from yeast (IsoDec, targets vs. decoys).
    /// These tests guard against circular self-validation by feeding the scorer
    /// envelopes whose peak intensities come from molecular-formula-based isotope
    /// distribution calculations, NOT from the Averagine model used internally.
    /// </summary>
    [TestFixture]
    public sealed class TestDeconvolutionScorerIndependentValidation
    {
        private static readonly AverageResidue Model = new Averagine();

        private static IsotopicEnvelope MakeEnvelope(
            double monoMass, int charge,
            IReadOnlyList<double> relativeIntensities,
            double baseIntensity = 1e6)
        {
            double monoMz = monoMass / charge + 1.00727646677;
            double isotopeStep = 1.003355 / charge;

            var peaks = new List<(double mz, double intensity)>();
            for (int i = 0; i < relativeIntensities.Count; i++)
            {
                double intensity = baseIntensity * relativeIntensities[i];
                if (intensity > 0)
                    peaks.Add((monoMz + i * isotopeStep, intensity));
            }

            return new IsotopicEnvelope(
                id: 0,
                peaks: peaks,
                monoisotopicmass: monoMass,
                chargestate: charge,
                intensity: peaks.Sum(p => p.intensity),
                score: 0.0);
        }

        // REGRESSION: failed before the fix — this test uses independently derived
        // isotope intensities that do NOT come from the Averagine model.
        [Test]
        public void ComputeFeatures_AngiotensinIIAtChargeTwo_CosineAboveThreshold()
        {
            // Angiotensin II: C₅₀H₇₁N₁₃O₁₂, mono mass 1045.5345, z=2
            // Formula-derived isotope distribution (normalized to M+0 = 1.0)
            var intensities = new[] { 1.0000, 0.5765, 0.1924, 0.0462, 0.0087 };
            var env = MakeEnvelope(1045.5345, 2, intensities);

            var features = DeconvolutionScorer.ComputeFeatures(env, Model);

            Assert.That(features.AveragineCosineSimilarity, Is.GreaterThan(0.85),
                $"Cosine should be > 0.85 for a real peptide envelope. Got {features.AveragineCosineSimilarity:F4}");
        }

        [Test]
        public void ComputeFeatures_AngiotensinIIAtChargeTwo_PpmErrorBelowFive()
        {
            var intensities = new[] { 1.0000, 0.5765, 0.1924, 0.0462, 0.0087 };
            var env = MakeEnvelope(1045.5345, 2, intensities);

            var features = DeconvolutionScorer.ComputeFeatures(env, Model);

            Assert.That(features.AvgPpmError, Is.LessThan(5.0),
                $"Exact m/z positions should yield ppm < 5. Got {features.AvgPpmError:F4}");
        }

        // REGRESSION: failed before the fix — validates alignment at a different mass range
        [Test]
        public void ComputeFeatures_InsulinBChainAtChargeFour_CosineAboveThreshold()
        {
            // Insulin B chain: C₁₅₇H₂₃₂N₄₂O₄₈S₁, mono mass 3495.6431, z=4
            // Formula-derived, normalized to most abundant (M+1) = 1.0
            var intensities = new[] { 0.5720, 1.0000, 0.9535, 0.6406, 0.3423, 0.1516, 0.0587 };
            var env = MakeEnvelope(3495.6431, 4, intensities);

            var features = DeconvolutionScorer.ComputeFeatures(env, Model);

            Assert.That(features.AveragineCosineSimilarity, Is.GreaterThan(0.85),
                $"Cosine should be > 0.85 for insulin B-chain. Got {features.AveragineCosineSimilarity:F4}");
        }

        [Test]
        public void ComputeFeatures_InsulinBChainAtChargeFour_CompletenessAboveHalf()
        {
            var intensities = new[] { 0.5720, 1.0000, 0.9535, 0.6406, 0.3423, 0.1516, 0.0587 };
            var env = MakeEnvelope(3495.6431, 4, intensities);

            var features = DeconvolutionScorer.ComputeFeatures(env, Model);

            Assert.That(features.PeakCompleteness, Is.GreaterThan(0.5),
                $"7 peaks should yield completeness > 0.5. Got {features.PeakCompleteness:F4}");
        }

        [Test]
        public void ComputeFeatures_IndependentEnvelope_OffByOneShift_CosineDrops()
        {
            // Edge case: simulate the exact off-by-one alignment bug the circular
            // tests cannot catch — shift all peaks by one isotope position.
            const double monoMass = 1045.5345;
            const int charge = 2;
            double monoMz = monoMass / charge + 1.00727646677;
            double isotopeStep = 1.003355 / charge;

            // Correct intensities but placed starting at M+1 position instead of M+0
            var intensities = new[] { 1.0000, 0.5765, 0.1924, 0.0462, 0.0087 };
            var peaks = new List<(double mz, double intensity)>();
            for (int i = 0; i < intensities.Length; i++)
            {
                peaks.Add((monoMz + (i + 1) * isotopeStep, 1e6 * intensities[i]));
            }

            var env = new IsotopicEnvelope(
                id: 0,
                peaks: peaks,
                monoisotopicmass: monoMass,
                chargestate: charge,
                intensity: peaks.Sum(p => p.intensity),
                score: 0.0);

            var features = DeconvolutionScorer.ComputeFeatures(env, Model);

            Assert.That(features.AveragineCosineSimilarity, Is.LessThan(0.7),
                $"Off-by-one shifted envelope should drop cosine well below the positive-case 0.85 baseline. Got {features.AveragineCosineSimilarity:F4}");
        }

        [Test]
        public void ComputeScore_IndependentAngiotensinEnvelope_ScoreAboveHalf()
        {
            var intensities = new[] { 1.0000, 0.5765, 0.1924, 0.0462, 0.0087 };
            var env = MakeEnvelope(1045.5345, 2, intensities);

            var features = DeconvolutionScorer.ComputeFeatures(env, Model);
            double score = DeconvolutionScorer.ComputeScore(features);

            Assert.That(score, Is.GreaterThan(0.5),
                $"A real peptide envelope should score > 0.5. Got {score:F4}");
            Assert.That(score, Is.InRange(0.0, 1.0),
                $"Score must be in [0,1]. Got {score:F4}");
        }

        [TestCase(800.0, 2, new[] { 1.0, 0.45, 0.12, 0.02 })]
        [TestCase(2000.0, 3, new[] { 0.75, 1.0, 0.70, 0.35, 0.13, 0.04 })]
        [TestCase(5000.0, 5, new[] { 0.30, 0.65, 1.0, 0.95, 0.70, 0.42, 0.21, 0.09 })]
        public void ComputeFeatures_VariousMassRanges_AllFeaturesFinite(
            double monoMass, int charge, double[] intensities)
        {
            var env = MakeEnvelope(monoMass, charge, intensities);

            var features = DeconvolutionScorer.ComputeFeatures(env, Model);

            Assert.That(double.IsFinite(features.AveragineCosineSimilarity), Is.True,
                "AveragineCosineSimilarity must be finite");
            Assert.That(double.IsFinite(features.AvgPpmError), Is.True,
                "AvgPpmError must be finite");
            Assert.That(double.IsFinite(features.PeakCompleteness), Is.True,
                "PeakCompleteness must be finite");
            Assert.That(double.IsFinite(features.IntensityRatioConsistency), Is.True,
                "IntensityRatioConsistency must be finite");
        }
    }
}
