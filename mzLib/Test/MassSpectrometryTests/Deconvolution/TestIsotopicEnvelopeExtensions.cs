using System;
using MassSpectrometry;
using NUnit.Framework;
using static Test.MassSpectrometryTests.Deconvolution.DeconvolutionTestHelpers;

namespace Test.Deconvolution
{
    /// <summary>
    /// Unit tests for <see cref="IsotopicEnvelopeExtensions"/> and the small read-only
    /// conveniences <see cref="IsotopicEnvelope.HasGenericScore"/> and
    /// <see cref="IsotopicEnvelope.GenericOrFallbackScore"/>.
    /// </summary>
    [TestFixture]
    public sealed class TestIsotopicEnvelopeExtensions
    {

        // ══════════════════════════════════════════════════════════════════════
        // GetOrComputeGenericScore: behaviour and contracts
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void GetOrComputeGenericScore_NotYetSet_ComputesAndStashes()
        {
            var env = BuildPerfectEnvelope();
            var parallelEnv = BuildPerfectEnvelope();
            double expected = DeconvolutionScorer.ScoreEnvelope(parallelEnv, Model);

            Assert.That(env.GenericScore, Is.Null,
                "Precondition: a freshly built envelope must not have a GenericScore");

            double returned = env.GetOrComputeGenericScore(Model);

            Assert.That(env.GenericScore.HasValue, Is.True,
                "After call, GenericScore must be populated on the envelope");
            Assert.That(returned, Is.EqualTo(env.GenericScore.Value),
                "Returned value must equal the value cached on the envelope");
            Assert.That(returned, Is.InRange(0.0, 1.0),
                $"Computed generic score must be in [0,1]. Got {returned:F4}");
            Assert.That(returned, Is.EqualTo(expected),
                "Cache-miss path must return DeconvolutionScorer.ScoreEnvelope(env, model) verbatim");
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

        /// <summary>
        /// Happy-path coverage for the <see cref="DeconvolutionParameters"/> overload
        /// (the common downstream case in MetaMorpheus). Verifies cache-miss behaviour:
        /// GenericScore is populated, the returned value is in [0, 1], and it equals
        /// what the underlying (envelope, AverageResidue) overload would produce —
        /// guarding against accidental delegation regressions.
        /// </summary>
        [Test]
        public void GetOrComputeGenericScore_ParametersOverload_NotYetSet_ComputesAndDelegatesToModelOverload()
        {
            var envViaParameters = BuildPerfectEnvelope();
            var envViaModel = BuildPerfectEnvelope();
            var parameters = new ClassicDeconvolutionParameters(
                minCharge: 1, maxCharge: 10, deconPpm: 4.0, intensityRatio: 3.0,
                averageResidueModel: Model);

            Assert.That(envViaParameters.GenericScore, Is.Null,
                "Precondition: a freshly built envelope must not have a GenericScore");

            double returnedFromParameters = envViaParameters.GetOrComputeGenericScore(parameters);
            double returnedFromModel = envViaModel.GetOrComputeGenericScore(Model);

            Assert.That(envViaParameters.GenericScore.HasValue, Is.True,
                "After call, GenericScore must be populated on the envelope");
            Assert.That(returnedFromParameters, Is.EqualTo(envViaParameters.GenericScore.Value),
                "Returned value must equal the value cached on the envelope");
            Assert.That(returnedFromParameters, Is.InRange(0.0, 1.0),
                $"Computed generic score must be in [0,1]. Got {returnedFromParameters:F4}");
            Assert.That(returnedFromParameters, Is.EqualTo(returnedFromModel),
                "Parameters overload must delegate to the AverageResidue overload and produce the same score");
        }

        /// <summary>
        /// Cache-hit coverage for the <see cref="DeconvolutionParameters"/> overload.
        /// A pre-set GenericScore that the underlying scorer would never produce
        /// must be returned verbatim, proving the parameters overload honours the
        /// cache rather than recomputing.
        /// </summary>
        [Test]
        public void GetOrComputeGenericScore_ParametersOverload_AlreadySet_ReturnsCachedWithoutRecompute()
        {
            var env = BuildPerfectEnvelope();
            env.SetGenericScore(0.42);
            var parameters = new ClassicDeconvolutionParameters(
                minCharge: 1, maxCharge: 10, deconPpm: 4.0, intensityRatio: 3.0,
                averageResidueModel: Model);

            double returned = env.GetOrComputeGenericScore(parameters);

            Assert.That(returned, Is.EqualTo(0.42),
                "Cached GenericScore must be returned verbatim without recomputation");
            Assert.That(env.GenericScore, Is.EqualTo(0.42),
                "Cache value on the envelope must remain unchanged");
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
