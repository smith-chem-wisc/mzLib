using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.MassSpectrometryTests.Deconvolution
{
    [TestFixture]
    public sealed class TestScoreEnvelopesNullHandling
    {
        private static readonly AverageResidue Model = new Averagine();

        private static IsotopicEnvelope MakeEnvelope(
            double monoMass, int charge, params (double mz, double intensity)[] peaks)
        {
            return new IsotopicEnvelope(
                id: 0,
                peaks: peaks.ToList(),
                monoisotopicmass: monoMass,
                chargestate: charge,
                intensity: peaks.Sum(p => p.intensity),
                score: 0.0);
        }

        // REGRESSION: failed before the fix — a single null in the batch
        // caused ArgumentNullException from ComputeFeatures, aborting the enumeration.
        [Test]
        public void ScoreEnvelopes_SequenceContainsNull_SkipsNullAndScoresRest()
        {
            double monoMass = 1000.0;
            int charge = 2;
            double mz = monoMass.ToMz(charge);
            var validEnvelope = MakeEnvelope(monoMass, charge, (mz, 1000.0));

            var envelopes = new List<IsotopicEnvelope> { validEnvelope, null, validEnvelope };

            var results = DeconvolutionScorer.ScoreEnvelopes(envelopes, Model).ToList();

            Assert.That(results, Has.Count.EqualTo(2));
            Assert.That(results[0].Envelope, Is.SameAs(validEnvelope));
            Assert.That(results[1].Envelope, Is.SameAs(validEnvelope));
            Assert.That(results[0].Score, Is.InRange(0.0, 1.0));
            Assert.That(results[1].Score, Is.InRange(0.0, 1.0));
        }

        [Test]
        public void ScoreEnvelopes_AllNulls_ReturnsEmpty()
        {
            var envelopes = new List<IsotopicEnvelope> { null, null, null };

            var results = DeconvolutionScorer.ScoreEnvelopes(envelopes, Model).ToList();

            Assert.That(results, Is.Empty);
        }

        [Test]
        public void ScoreEnvelopes_NullCollection_ThrowsArgumentNullException()
        {
            // No .ToList() inside the lambda -- this only passes if the throw fires
            // at the call site, before MoveNext. Pins the wrapper+iterator split
            // against a silent regression to a single iterator method (which would
            // make the throw deferred to first enumeration).
            Assert.That(
                () => DeconvolutionScorer.ScoreEnvelopes(null, Model),
                Throws.TypeOf<ArgumentNullException>()
                    .With.Property("ParamName").EqualTo("envelopes"));
        }

        [Test]
        public void ScoreEnvelopes_NullModel_ThrowsArgumentNullException()
        {
            var envelopes = new List<IsotopicEnvelope>();

            // See comment on the null-collection sibling -- no .ToList() so this
            // distinguishes eager throw from deferred-to-MoveNext throw.
            Assert.That(
                () => DeconvolutionScorer.ScoreEnvelopes(envelopes, null),
                Throws.TypeOf<ArgumentNullException>()
                    .With.Property("ParamName").EqualTo("model"));
        }

        [Test]
        public void ScoreEnvelopes_EmptySequence_ReturnsEmpty()
        {
            var envelopes = Enumerable.Empty<IsotopicEnvelope>();

            var results = DeconvolutionScorer.ScoreEnvelopes(envelopes, Model).ToList();

            Assert.That(results, Is.Empty);
        }

        [Test]
        public void ScoreEnvelopes_ValidEnvelopes_ReturnsPairPerEnvelope()
        {
            double monoMass1 = 800.0;
            double monoMass2 = 1200.0;
            int charge = 3;

            var env1 = MakeEnvelope(monoMass1, charge, (monoMass1.ToMz(charge), 500.0));
            var env2 = MakeEnvelope(monoMass2, charge, (monoMass2.ToMz(charge), 750.0));

            var results = DeconvolutionScorer.ScoreEnvelopes(new[] { env1, env2 }, Model).ToList();

            Assert.That(results, Has.Count.EqualTo(2));
            Assert.That(results[0].Envelope, Is.SameAs(env1));
            Assert.That(results[1].Envelope, Is.SameAs(env2));
            Assert.That(results[0].Score, Is.InRange(0.0, 1.0));
            Assert.That(results[1].Score, Is.InRange(0.0, 1.0));
        }

        [Test]
        public void ScoreEnvelopes_NullAtStartAndEnd_SkipsAllNulls()
        {
            double monoMass = 950.0;
            int charge = 2;
            var env = MakeEnvelope(monoMass, charge, (monoMass.ToMz(charge), 600.0));

            var envelopes = new List<IsotopicEnvelope> { null, env, null };

            var results = DeconvolutionScorer.ScoreEnvelopes(envelopes, Model).ToList();

            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].Envelope, Is.SameAs(env));
        }
    }
}
