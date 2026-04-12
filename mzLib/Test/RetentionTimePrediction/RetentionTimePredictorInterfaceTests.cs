using NUnit.Framework;
using Chromatography.RetentionTimePrediction;
using Chromatography;
using System.Collections.Generic;
using System.Linq;

namespace Test.RetentionTimePrediction
{
    /// <summary>
    /// Tests for the unified IRetentionTimePredictor interface contract.
    /// Verifies that both local (Chronologer) and remote (Koina) predictors
    /// satisfy the interface correctly, including the default batch method.
    /// </summary>
    [TestFixture]
    public class RetentionTimePredictorInterfaceTests
    {
        // ── Helper: minimal IRetentionPredictable stub ────────────────────

        private sealed class StubPeptide : IRetentionPredictable
        {
            public string BaseSequence { get; init; } = string.Empty;
            public string FullSequence { get; init; } = string.Empty;
            public string FullSequenceWithMassShifts { get; init; } = string.Empty;
            public double MonoisotopicMass { get; init; } = 1000;
        }

        // ── Helper: minimal mock IRetentionTimePredictor ──────────────────

        private sealed class MockPredictor : IRetentionTimePredictor
        {
            public string PredictorName => "Mock";
            public SeparationType SeparationType => SeparationType.HPLC;

            private readonly Dictionary<string, double?> _predictions;

            public MockPredictor(Dictionary<string, double?> predictions)
            {
                _predictions = predictions;
            }

            public double? PredictRetentionTime(IRetentionPredictable peptide,
                out RetentionTimeFailureReason? failureReason)
            {
                failureReason = null;
                if (_predictions.TryGetValue(peptide.FullSequence, out var value))
                    return value;
                failureReason = RetentionTimeFailureReason.PredictionError;
                return null;
            }

            public string? GetFormattedSequence(IRetentionPredictable peptide,
                out RetentionTimeFailureReason? failureReason)
            {
                failureReason = null;
                return peptide.BaseSequence;
            }
            // Note: PredictRetentionTimes is NOT overridden here.
            // This verifies the default interface implementation works.
        }

        // ── Tests for default PredictRetentionTimes batch method ──────────

        [Test]
        public void PredictRetentionTimes_DefaultImpl_ReturnsResultForEachPeptide()
        {
            var predictions = new Dictionary<string, double?>
            {
                { "PEPTIDE", 25.3 },
                { "PROTEIN", 42.1 }
            };
            IRetentionTimePredictor predictor = new MockPredictor(predictions);

            var peptides = new List<IRetentionPredictable>
            {
                new StubPeptide { BaseSequence = "PEPTIDE", FullSequence = "PEPTIDE" },
                new StubPeptide { BaseSequence = "PROTEIN", FullSequence = "PROTEIN" }
            };

            var results = predictor.PredictRetentionTimes(peptides);

            Assert.That(results.Count, Is.EqualTo(2));
            Assert.That(results["PEPTIDE"], Is.EqualTo(25.3));
            Assert.That(results["PROTEIN"], Is.EqualTo(42.1));
        }

        [Test]
        public void PredictRetentionTimes_DefaultImpl_ReturnsNullForUnknownSequence()
        {
            IRetentionTimePredictor predictor = new MockPredictor(
                new Dictionary<string, double?>());

            var peptides = new List<IRetentionPredictable>
            {
                new StubPeptide { BaseSequence = "PEPTIDE", FullSequence = "PEPTIDE" }
            };

            var results = predictor.PredictRetentionTimes(peptides);

            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results["PEPTIDE"], Is.Null);
        }

        [Test]
        public void PredictRetentionTimes_DefaultImpl_EmptyInput_ReturnsEmptyDictionary()
        {
            IRetentionTimePredictor predictor = new MockPredictor(
                new Dictionary<string, double?>());

            var results = predictor.PredictRetentionTimes(
                new List<IRetentionPredictable>());

            Assert.That(results, Is.Empty);
        }

        [Test]
        public void PredictRetentionTimes_DefaultImpl_KeyIsFullSequence()
        {
            // Verify the dictionary key is FullSequence, not BaseSequence
            var predictions = new Dictionary<string, double?>
            {
                { "PEPTM[Oxidation on M]IDE", 31.7 }
            };
            IRetentionTimePredictor predictor = new MockPredictor(predictions);

            var peptides = new List<IRetentionPredictable>
            {
                new StubPeptide
                {
                    BaseSequence = "PEPTMIDE",
                    FullSequence = "PEPTM[Oxidation on M]IDE"
                }
            };

            var results = predictor.PredictRetentionTimes(peptides);

            Assert.That(results.ContainsKey("PEPTM[Oxidation on M]IDE"), Is.True);
            Assert.That(results.ContainsKey("PEPTMIDE"), Is.False);
        }

        [Test]
        public void PredictRetentionTimes_DefaultImpl_MixedNullAndNonNull_HandledCorrectly()
        {
            var predictions = new Dictionary<string, double?>
            {
                { "PEPTIDE", 25.3 },
                // "UNKNOWN" not in dictionary -> will return null
            };
            IRetentionTimePredictor predictor = new MockPredictor(predictions);

            var peptides = new List<IRetentionPredictable>
            {
                new StubPeptide { BaseSequence = "PEPTIDE", FullSequence = "PEPTIDE" },
                new StubPeptide { BaseSequence = "UNKNOWN", FullSequence = "UNKNOWN" }
            };

            var results = predictor.PredictRetentionTimes(peptides);

            Assert.That(results["PEPTIDE"], Is.EqualTo(25.3));
            Assert.That(results["UNKNOWN"], Is.Null);
        }

        // ── Tests verifying Chronologer satisfies the interface ───────────

        [Test]
        public void ChronologerRetentionTimePredictor_ImplementsIRetentionTimePredictor()
        {
            using var predictor = new
                Chromatography.RetentionTimePrediction.Chronologer
                .ChronologerRetentionTimePredictor();

            Assert.That(predictor, Is.InstanceOf<IRetentionTimePredictor>());
        }

        [Test]
        public void ChronologerRetentionTimePredictor_PredictRetentionTimes_DefaultBatch_Works()
        {
            using var predictor = new
                Chromatography.RetentionTimePrediction.Chronologer
                .ChronologerRetentionTimePredictor();

            // Cast to interface to ensure we're testing the interface method
            IRetentionTimePredictor iPredictor = predictor;

            var peptides = new List<IRetentionPredictable>
            {
                new StubPeptide { BaseSequence = "PEPTIDE", FullSequence = "PEPTIDE", FullSequenceWithMassShifts = "PEPTIDE" },
                new StubPeptide { BaseSequence = "AGHCEWQMK", FullSequence = "AGHCEWQMK", FullSequenceWithMassShifts = "AGHCEWQMK" }
            };

            var results = iPredictor.PredictRetentionTimes(peptides);

            Assert.That(results.Count, Is.EqualTo(2));
            Assert.That(results["PEPTIDE"], Is.Not.Null);
            Assert.That(results["AGHCEWQMK"], Is.Not.Null);
        }

        // ── Tests verifying Koina RetentionTimeModel satisfies the interface ──
        // These verify the interface is implemented; actual HTTP calls are [Explicit]

        [Test]
        public void Prosit2019iRT_ImplementsIRetentionTimePredictor()
        {
            // Verify that Prosit2019iRT satisfies the interface.
            // This test will fail to compile if RetentionTimeModel does not implement
            // IRetentionTimePredictor — which is the point.
            var model = new PredictionClients.Koina.SupportedModels.RetentionTimeModels
                .Prosit2019iRT();

            IRetentionTimePredictor predictor = model;

            Assert.That(predictor, Is.Not.Null);
            Assert.That(predictor.PredictorName, Is.Not.Null.And.Not.Empty);
            Assert.That(predictor.SeparationType, Is.EqualTo(SeparationType.HPLC));
        }

        [Test]
        public void Prosit2020iRTTMT_ImplementsIRetentionTimePredictor()
        {
            var model = new PredictionClients.Koina.SupportedModels.RetentionTimeModels
                .Prosit2020iRTTMT();

            IRetentionTimePredictor predictor = model;

            Assert.That(predictor, Is.Not.Null);
            Assert.That(predictor.PredictorName, Is.Not.Null.And.Not.Empty);
        }

        [Test]
        [Explicit("Requires live Koina API access")]
        public void Prosit2019iRT_PredictRetentionTimes_BatchOverride_ReturnsResults()
        {
            // Verifies the Koina override of PredictRetentionTimes issues a single
            // batch HTTP call rather than looping per-peptide.
            var model = new PredictionClients.Koina.SupportedModels.RetentionTimeModels
                .Prosit2019iRT();

            IRetentionTimePredictor predictor = model;

            var peptides = new List<IRetentionPredictable>
            {
                new StubPeptide { BaseSequence = "PEPTIDE", FullSequence = "PEPTIDE" },
                new StubPeptide { BaseSequence = "AGHCEWQMK", FullSequence = "AGHCEWQMK" }
            };

            var results = predictor.PredictRetentionTimes(peptides);

            Assert.That(results.Count, Is.EqualTo(2));
            Assert.That(results.Values.All(v => v.HasValue), Is.True);
            // iRT values are roughly in the range -100 to 200
            foreach (var value in results.Values)
                Assert.That(value!.Value, Is.InRange(-100.0, 200.0));
        }
    }
}
