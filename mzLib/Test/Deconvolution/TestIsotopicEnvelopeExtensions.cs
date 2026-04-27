using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.Deconvolution
{
    /// <summary>
    /// Unit tests for <see cref="IsotopicEnvelopeExtensions"/> and the small read-only
    /// conveniences <see cref="IsotopicEnvelope.HasGenericScore"/> and
    /// <see cref="IsotopicEnvelope.GenericOrFallbackScore"/>. Mirrors the envelope-construction
    /// pattern from <see cref="TestDeconvolutionScorerUnit"/> (6-arg constructor with named args)
    /// to avoid introducing parallel helpers.
    /// </summary>
    [TestFixture]
    public sealed class TestIsotopicEnvelopeExtensions
    {
        // ── Shared model and test mass ────────────────────────────────────────

        private static readonly AverageResidue Model = new Averagine();

        /// <summary>
        /// ~5 kDa peptide — well within the Averagine model's calibrated range,
        /// produces a clearly resolved isotope pattern at z = 5. Matches
        /// <see cref="TestDeconvolutionScorerUnit"/>.
        /// </summary>
        private const double TestMass = 5000.0;
        private const int TestCharge = 5;

        // ── Helper (mirrors TestDeconvolutionScorerUnit.BuildPerfectEnvelope) ─

        /// <summary>
        /// Builds a perfect synthetic envelope using the same 6-arg constructor pattern with
        /// named args as <see cref="TestDeconvolutionScorerUnit"/>. The optional
        /// <paramref name="score"/> parameter lets tests inject a known algorithm-specific
        /// score for fallback verification.
        /// </summary>
        private static IsotopicEnvelope BuildPerfectEnvelope(
            double monoMass = TestMass,
            int charge = TestCharge,
            double baseIntens = 1e6,
            double score = 0.999)
        {
            int avgIdx = Model.GetMostIntenseMassIndex(monoMass);

            double[] rawMasses = Model.GetAllTheoreticalMasses(avgIdx);
            double[] rawIntens = Model.GetAllTheoreticalIntensities(avgIdx);

            // Mass-ascending sort so index 0 = monoisotopic
            var sorted = rawMasses.Zip(rawIntens)
                .OrderBy(p => p.First)
                .ToArray();

            int absCharge = Math.Abs(charge);
            double isotopeStep = Constants.C13MinusC12 / absCharge;
            double monoMz = monoMass.ToMz(charge);
            var peaks = new List<(double mz, double intensity)>();

            for (int n = 0; n < sorted.Length; n++)
            {
                double intensity = baseIntens * sorted[n].Second;
                if (intensity < baseIntens * 0.001) continue; // skip near-zero peaks
                double mz = monoMz + n * isotopeStep;
                peaks.Add((mz, intensity));
            }

            return new IsotopicEnvelope(
                id: 0,
                peaks: peaks,
                monoisotopicmass: monoMass,
                chargestate: charge,
                intensity: peaks.Sum(p => p.intensity),
                score: score);
        }

        // ══════════════════════════════════════════════════════════════════════
        // GetOrComputeGenericScore: behaviour and contracts
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void GetOrComputeGenericScore_NotYetSet_ComputesAndStashes()
        {
            var env = BuildPerfectEnvelope();
            Assert.That(env.GenericScore, Is.Null,
                "Precondition: a freshly built envelope must not have a GenericScore");

            double returned = env.GetOrComputeGenericScore(Model);

            Assert.That(env.GenericScore.HasValue, Is.True,
                "After call, GenericScore must be populated on the envelope");
            Assert.That(returned, Is.EqualTo(env.GenericScore.Value),
                "Returned value must equal the value cached on the envelope");
            Assert.That(returned, Is.InRange(0.0, 1.0),
                $"Computed generic score must be in [0,1]. Got {returned:F4}");
        }

        [Test]
        public void GetOrComputeGenericScore_AlreadySet_ReturnsCachedWithoutRecompute()
        {
            // Build an envelope that would NOT score 0.42 if recomputed — that way, getting
            // 0.42 back is proof the cache was honoured rather than overwritten.
            var env = BuildPerfectEnvelope();
            env.SetGenericScore(0.42);

            double returned = env.GetOrComputeGenericScore(Model);

            Assert.That(returned, Is.EqualTo(0.42),
                "Cached GenericScore must be returned verbatim without recomputation");
            Assert.That(env.GenericScore, Is.EqualTo(0.42),
                "Cache value on the envelope must remain unchanged");
        }

        [Test]
        public void GetOrComputeGenericScore_NullEnvelope_Throws()
        {
            IsotopicEnvelope nullEnv = null;

            var ex = Assert.Throws<ArgumentNullException>(
                () => nullEnv.GetOrComputeGenericScore(Model));
            Assert.That(ex.ParamName, Is.EqualTo("envelope"));
        }

        [Test]
        public void GetOrComputeGenericScore_NullModel_Throws()
        {
            var env = BuildPerfectEnvelope();

            // (envelope, AverageResidue) overload with null model
            var ex = Assert.Throws<ArgumentNullException>(
                () => env.GetOrComputeGenericScore((AverageResidue)null));
            Assert.That(ex.ParamName, Is.EqualTo("model"));

            // (envelope, DeconvolutionParameters) overload with null parameters
            var ex2 = Assert.Throws<ArgumentNullException>(
                () => env.GetOrComputeGenericScore((DeconvolutionParameters)null));
            Assert.That(ex2.ParamName, Is.EqualTo("parameters"));
        }

        // ══════════════════════════════════════════════════════════════════════
        // HasGenericScore / GenericOrFallbackScore: round-trip and fallback
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void HasGenericScore_TrueAfterSet()
        {
            var env = BuildPerfectEnvelope();

            Assert.That(env.HasGenericScore, Is.False,
                "HasGenericScore must be false before any score is set");

            env.SetGenericScore(0.75);

            Assert.That(env.HasGenericScore, Is.True,
                "HasGenericScore must be true after SetGenericScore");
        }

        [Test]
        public void GenericOrFallbackScore_FallsBackToScoreWhenNotSet()
        {
            const double knownScore = 0.31415;
            var env = BuildPerfectEnvelope(score: knownScore);

            Assert.That(env.GenericScore, Is.Null,
                "Precondition: GenericScore must be unset");
            Assert.That(env.GenericOrFallbackScore, Is.EqualTo(knownScore),
                "When GenericScore is unset, GenericOrFallbackScore must return the algorithm-specific Score");
        }

        [Test]
        public void GenericOrFallbackScore_ReturnsGenericWhenSet()
        {
            const double algoScore = 0.123;
            const double genericScore = 0.876;

            var env = BuildPerfectEnvelope(score: algoScore);
            env.SetGenericScore(genericScore);

            Assert.That(env.GenericOrFallbackScore, Is.EqualTo(genericScore),
                "When GenericScore is set, GenericOrFallbackScore must return it (not Score)");
            Assert.That(env.Score, Is.EqualTo(algoScore),
                "Score itself must remain the algorithm-specific value");
        }
    }
}
