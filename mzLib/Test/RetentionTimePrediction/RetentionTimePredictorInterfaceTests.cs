using NUnit.Framework;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.Chronologer;
using Chromatography;
using System.Collections.Generic;
using System.Linq;
using Omics.SequenceConversion;

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
        // ── Shared constants ──────────────────────────────────────────────

        private const double StubMonoisotopicMass = 1000.0;

        // ── Helper: minimal IRetentionPredictable stub ────────────────────

        private sealed class StubPeptide : IRetentionPredictable
        {
            public string BaseSequence { get; init; } = string.Empty;
            public string FullSequence { get; init; } = string.Empty;
            public string FullSequenceWithMassShifts { get; init; } = string.Empty;
            public double MonoisotopicMass { get; init; } = StubMonoisotopicMass;
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

            // For the mock, equivalent == time. Delegate to PredictRetentionTime.
            public double? PredictRetentionTimeEquivalent(IRetentionPredictable peptide,
                out RetentionTimeFailureReason? failureReason)
                => PredictRetentionTime(peptide, out failureReason);

            // Minimal batch stub: per-peptide via the single-shot path.
            public IReadOnlyList<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)>
                PredictRetentionTimeEquivalents(IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1)
            {
                var output = new List<(double?, IRetentionPredictable, RetentionTimeFailureReason?)>();
                foreach (var p in peptides)
                {
                    if (p is null) continue;
                    var v = PredictRetentionTime(p, out var reason);
                    output.Add((v, p, reason));
                }
                return output;
            }

            public void Dispose() { }

            // PredictRetentionTimes is NOT overridden here —
            // verifies the default interface implementation works correctly.
        }

        // ── Chronologer shared instance (expensive to construct, lazy) ────
        //
        // Lazy init isolates TorchSharp model load to the three Chronologer-using
        // tests; failures there do not collateral-fail the mock-only tests in
        // this fixture.

        private static ChronologerRetentionTimePredictor? _chronologer;

        private static ChronologerRetentionTimePredictor Chronologer =>
            _chronologer ??= new ChronologerRetentionTimePredictor(
                SequenceConversionHandlingMode.RemoveIncompatibleElements);

        [OneTimeTearDown]
        public static void OneTimeTearDown()
        {
            _chronologer?.Dispose();
        }

        // ── Tests: default PredictRetentionTimes batch method ─────────────

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
        public void PredictRetentionTimes_DefaultImpl_NullInput_ThrowsArgumentNullException()
        {
            IRetentionTimePredictor predictor = new MockPredictor(
                new Dictionary<string, double?>());

            Assert.That(() => predictor.PredictRetentionTimes(null!),
                Throws.ArgumentNullException);
        }

        [Test]
        public void PredictRetentionTimes_DefaultImpl_KeyIsFullSequence()
        {
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
                { "PEPTIDE", 25.3 }
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

        [Test]
        public void PredictRetentionTimes_DefaultImpl_ReturnsIReadOnlyDictionary()
        {
            IRetentionTimePredictor predictor = new MockPredictor(
                new Dictionary<string, double?> { { "PEPTIDE", 1.0 } });

            var results = predictor.PredictRetentionTimes(new List<IRetentionPredictable>
            {
                new StubPeptide { FullSequence = "PEPTIDE" }
            });

            Assert.That(results, Is.InstanceOf<IReadOnlyDictionary<string, double?>>());
            Assert.That(results as Dictionary<string, double?>, Is.Null,
                "Result should not be a mutable Dictionary");
        }

        // ── Tests: failureReason is not propagated in batch (documented limitation) ──

        [Test]
        public void PredictRetentionTimes_DefaultImpl_NullResultEntries_FailureReasonNotAvailable()
        {
            IRetentionTimePredictor predictor = new MockPredictor(
                new Dictionary<string, double?>());

            var results = predictor.PredictRetentionTimes(new List<IRetentionPredictable>
            {
                new StubPeptide { FullSequence = "UNKNOWN" }
            });

            Assert.That(results["UNKNOWN"], Is.Null);
        }

        // ── Tests: Chronologer satisfies the interface ────────────────────

        [Test]
        public void ChronologerRetentionTimePredictor_ImplementsIRetentionTimePredictor()
        {
            Assert.That(Chronologer, Is.InstanceOf<IRetentionTimePredictor>());
        }

        [Test]
        public void ChronologerRetentionTimePredictor_PredictRetentionTimes_DefaultBatch_Works()
        {
            IRetentionTimePredictor iPredictor = Chronologer;

            var peptides = new List<IRetentionPredictable>
            {
                new StubPeptide { BaseSequence = "PEPTIDE", FullSequence = "PEPTIDE" },
                new StubPeptide { BaseSequence = "AGHCEWQMK", FullSequence = "AGHCEWQMK" }
            };

            var results = iPredictor.PredictRetentionTimes(peptides);

            Assert.That(results.Count, Is.EqualTo(2));
            Assert.That(results["PEPTIDE"], Is.Not.Null);
            Assert.That(results["AGHCEWQMK"], Is.Not.Null);
        }

        [Test]
        public void ChronologerRetentionTimePredictor_PredictRetentionTimes_ReturnsIReadOnlyDictionary()
        {
            IRetentionTimePredictor iPredictor = Chronologer;

            var results = iPredictor.PredictRetentionTimes(new List<IRetentionPredictable>
            {
                new StubPeptide { BaseSequence = "PEPTIDE", FullSequence = "PEPTIDE" }
            });

            Assert.That(results, Is.InstanceOf<IReadOnlyDictionary<string, double?>>());
        }

        // ── Tests: Koina models satisfy the interface ─────────────────────

        [Test]
        public void Prosit2019iRT_ImplementsIRetentionTimePredictor()
        {
            IRetentionTimePredictor predictor = new
                PredictionClients.Koina.SupportedModels.RetentionTimeModels.Prosit2019iRT();

            Assert.That(predictor, Is.Not.Null);
            Assert.That(predictor.PredictorName, Is.Not.Null.And.Not.Empty);
            Assert.That(predictor.SeparationType, Is.EqualTo(SeparationType.HPLC));
        }

        [Test]
        public void Prosit2020iRTTMT_ImplementsIRetentionTimePredictor()
        {
            IRetentionTimePredictor predictor = new
                PredictionClients.Koina.SupportedModels.RetentionTimeModels
                .Prosit2020iRTTMT();

            Assert.That(predictor, Is.Not.Null);
            Assert.That(predictor.PredictorName, Is.Not.Null.And.Not.Empty);
        }

        [Test]
        [Category("Integration")]
        public void Prosit2019iRT_PredictRetentionTime_SinglePeptide_DoesNotThrow()
        {
            IRetentionTimePredictor predictor = new
                PredictionClients.Koina.SupportedModels.RetentionTimeModels.Prosit2019iRT();

            Assert.DoesNotThrow(() =>
            {
                predictor.PredictRetentionTime(
                    new StubPeptide { BaseSequence = "PEPTIDE", FullSequence = "PEPTIDE" },
                    out _);
            });
        }

        [Test]
        public void Prosit2019iRT_GetFormattedSequence_ReturnsNonNullForValidPeptide()
        {
            IRetentionTimePredictor predictor = new
                PredictionClients.Koina.SupportedModels.RetentionTimeModels.Prosit2019iRT();

            var result = predictor.GetFormattedSequence(
                new StubPeptide { BaseSequence = "PEPTIDE", FullSequence = "PEPTIDE" },
                out var failureReason);

            Assert.That(result, Is.Not.Null);
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void Prosit2019iRT_GetFormattedSequence_NullPeptide_ReturnsNullWithFailureReason()
        {
            IRetentionTimePredictor predictor = new
                PredictionClients.Koina.SupportedModels.RetentionTimeModels.Prosit2019iRT();

            var result = predictor.GetFormattedSequence(null!, out var failureReason);

            Assert.That(result, Is.Null);
            Assert.That(failureReason, Is.Not.Null);
        }

        // ── Live API tests (run in CI under the Integration category) ────

        [Test]
        [Category("Integration")]
        public void Prosit2019iRT_PredictRetentionTimes_BatchOverride_ReturnsResults()
        {
            IRetentionTimePredictor predictor = new
                PredictionClients.Koina.SupportedModels.RetentionTimeModels.Prosit2019iRT();

            var peptides = new List<IRetentionPredictable>
            {
                new StubPeptide { BaseSequence = "PEPTIDE", FullSequence = "PEPTIDE" },
                new StubPeptide { BaseSequence = "AGHCEWQMK", FullSequence = "AGHCEWQMK" }
            };

            var results = predictor.PredictRetentionTimes(peptides);

            Assert.That(results.Count, Is.EqualTo(2));
            Assert.That(results.Values.All(v => v.HasValue), Is.True);
            foreach (var value in results.Values)
                Assert.That(value!.Value, Is.InRange(-100.0, 200.0));
        }

        [Test]
        [Category("Integration")]
        public void Prosit2019iRT_PredictRetentionTimes_MixedValidInvalid_InvalidReturnsNull()
        {
            // Regression test for the pre-existing ModelInputs index bug.
            IRetentionTimePredictor predictor = new
                PredictionClients.Koina.SupportedModels.RetentionTimeModels.Prosit2019iRT();

            var peptides = new List<IRetentionPredictable>
            {
                new StubPeptide { BaseSequence = "PEPUIDE", FullSequence = "PEPUIDE" },
                new StubPeptide { BaseSequence = "AGHCEWQMK", FullSequence = "AGHCEWQMK" }
            };

            var results = predictor.PredictRetentionTimes(peptides);

            Assert.That(results["PEPUIDE"], Is.Null,
                "Invalid sequence should produce null prediction");
            Assert.That(results["AGHCEWQMK"], Is.Not.Null,
                "Valid sequence after invalid one must receive correct prediction");
        }
    }
}
