using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Threading.Tasks;
using Chromatography;
using Chromatography.RetentionTimePrediction;
using MzLibUtil;
using NUnit.Framework;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;

namespace Test.KoinaTests.RetentionTimePrediction
{
    /// <summary>
    /// Unit tests for <see cref="RetentionTimeModel"/> that run without any network access.
    ///
    /// Coverage strategy:
    ///   • A test double (<see cref="TestableRetentionTimeModel"/>) subclasses
    ///     <see cref="RetentionTimeModel"/> and overrides <c>AsyncThrottledPredictor</c>
    ///     so every code path that would normally call the Koina HTTP endpoint can be
    ///     exercised with injected fake predictions.
    ///   • Paths that sit before/after the HTTP call (input validation, deduplication,
    ///     result realignment, JSON response parsing) can therefore be tested directly.
    ///   • The HTTP section itself (~30 lines of <c>AsyncThrottledPredictor</c>) is not
    ///     unit-testable without a refactor; it is covered by the Integration-category
    ///     tests that hit the real Koina API.
    ///   • Three IRetentionTimePredictor members are intentionally minimal and not
    ///     unit-tested directly: <c>PredictRetentionTimeEquivalent</c> (a one-line
    ///     delegation to <c>PredictRetentionTime</c> with no behavior of its own),
    ///     <c>Dispose</c> (empty no-op default; subclasses with unmanaged resources
    ///     override), and the multi-batch throttling delay inside
    ///     <c>AsyncThrottledPredictor</c> (a single <c>await Task.Delay</c> that only
    ///     fires when an input batch spans more than one chunk; pinning it would
    ///     require timing-aware assertions that don't add behavioral value).
    ///
    /// These tests intentionally live OUTSIDE the [Category("Integration")] fixtures so
    /// they run on every PR build.
    /// </summary>
    [TestFixture]
    public class RetentionTimeModelTests
    {
        // ── Shared Unimod-35 / Unimod-4 converter (same as Prosit2019iRT) ─────────
        //
        // Using a real converter lets TryCleanSequence exercise its full parse /
        // serialize path for canonical and non-canonical sequences.

        private static readonly IReadOnlySet<int> SupportedMods = new HashSet<int> { 35, 4 };

        // ── Minimal IRetentionPredictable stub ────────────────────────────────────

        private sealed class StubPeptide : IRetentionPredictable
        {
            public string BaseSequence { get; init; } = string.Empty;
            public string FullSequence { get; init; } = string.Empty;
            public string FullSequenceWithMassShifts { get; init; } = string.Empty;
            public double MonoisotopicMass { get; init; } = 1000.0;
        }

        // ── Test double: bypasses HTTP, exposes protected hooks ───────────────────
        //
        // When BypassHttpWithStub is true, AsyncThrottledPredictor invokes the
        // configured PredictStub (or throws, if ThrowOnPredict is set) instead of
        // calling the real base implementation. Set BypassHttpWithStub to false to
        // run the real base logic — safe only when modelInputs is null/empty, or
        // when every input is guaranteed to fail validation so HTTP is skipped.

        private sealed class TestableRetentionTimeModel : RetentionTimeModel
        {
            public TestableRetentionTimeModel(
                SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.RemoveIncompatibleElements)
                : base(CreateUnimodConverter(UnimodSequenceFormatSchema.Instance, SupportedMods))
            {
                ModHandlingMode = mode;
                MaxNumberOfBatchesPerRequest = 500;
                ThrottlingDelayInMilliseconds = 0;
            }

            public override string ModelName => "Test_RT_Model";
            public override int MaxBatchSize => 1000;
            public override int MaxNumberOfBatchesPerRequest { get; init; }
            public override int ThrottlingDelayInMilliseconds { get; init; }
            public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 100;
            public override SequenceConversionHandlingMode ModHandlingMode { get; init; }
            public override int MaxPeptideLength => 30;
            public override int MinPeptideLength => 1;
            public override bool IsIndexedRetentionTimeModel => true;
            public override IReadOnlySet<int> AllowedUnimodIds => SupportedMods;

            // Test knobs
            public bool BypassHttpWithStub { get; set; } = true;
            public bool ThrowOnPredict { get; set; }
            public Func<List<RetentionTimePredictionInput>, List<PeptideRTPrediction>>? PredictStub { get; set; }

            protected override List<Dictionary<string, object>> ToBatchedRequests(
                List<RetentionTimePredictionInput> validInputs)
                => new();

            protected override Task<List<PeptideRTPrediction>> AsyncThrottledPredictor(
                List<RetentionTimePredictionInput> modelInputs)
            {
                if (!BypassHttpWithStub)
                    return base.AsyncThrottledPredictor(modelInputs);

                if (ThrowOnPredict)
                    throw new Exception("Simulated prediction failure");

                var stubResults = PredictStub?.Invoke(modelInputs ?? new List<RetentionTimePredictionInput>())
                                  ?? new List<PeptideRTPrediction>();
                Predictions = stubResults;
                return Task.FromResult(stubResults);
            }

            // Public proxy so tests can hit the protected method directly without
            // having to round-trip through Predict + a full response pipeline.
            public List<PeptideRTPrediction> InvokeResponseToPredictions(
                IReadOnlyList<string> responses,
                List<RetentionTimePredictionInput> requestInputs)
                => ResponseToPredictions(responses, requestInputs);
        }

        // ── Helpers ───────────────────────────────────────────────────────────────

        private static PeptideRTPrediction Pred(string fullSeq, double? rt, bool? indexed = true, string? validated = null)
            => new(fullSeq, validated ?? fullSeq, rt, indexed);

        private static string BuildJsonResponse(params double[] rtValues)
        {
            // Minimal ResponseJSONStruct-compatible JSON: single output with a Data list.
            var dataList = string.Join(",", rtValues.Select(v => v.ToString(System.Globalization.CultureInfo.InvariantCulture)));
            return $"{{\"outputs\":[{{\"name\":\"rt\",\"datatype\":\"FP32\",\"shape\":[{rtValues.Length}],\"data\":[{dataList}]}}]}}";
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // IRetentionTimePredictor properties
        // ═══════════════════════════════════════════════════════════════════════════

        /// <summary>
        /// PredictorName must surface the Koina ModelName so callers can tell which
        /// model produced a given prediction. This is the property MetaMorpheus uses
        /// to match PSMs against `CommonParameters.RTPredictorName`.
        /// </summary>
        [Test]
        public void PredictorName_DelegatesToModelName()
        {
            var model = new TestableRetentionTimeModel();
            Assert.That(model.PredictorName, Is.EqualTo("Test_RT_Model"));
            Assert.That(model.PredictorName, Is.EqualTo(model.ModelName));
        }

        /// <summary>
        /// All Koina RT models target HPLC by default. The property is virtual so a
        /// future CE-specific Koina model could override it; this test pins the
        /// default so any accidental change becomes visible.
        /// </summary>
        [Test]
        public void SeparationType_DefaultsToHPLC()
        {
            var model = new TestableRetentionTimeModel();
            Assert.That(model.SeparationType, Is.EqualTo(SeparationType.HPLC));
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // PredictRetentionTime (single-peptide, delegates to batch via list-of-one)
        // ═══════════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Null peptide must not NRE -- it must return null with EmptySequence as the
        /// failure reason. A null peptide is semantically "no input to predict from"
        /// (same bucket as an empty BaseSequence/FullSequence), not "the model was
        /// consulted and failed". EmptySequence lets upstream callers route null-input
        /// cases distinctly from genuine model errors. PEP downstream treats null as a
        /// sentinel value; propagating an NRE would crash the whole analysis.
        /// </summary>
        [Test]
        public void PredictRetentionTime_NullPeptide_ReturnsNullWithEmptySequence()
        {
            var model = new TestableRetentionTimeModel();

            var result = model.PredictRetentionTime(null!, out var reason);

            Assert.That(result, Is.Null);
            Assert.That(reason, Is.EqualTo(RetentionTimeFailureReason.EmptySequence));
        }

        /// <summary>
        /// An empty BaseSequence is semantically different from a prediction failure:
        /// there is nothing to predict at all. The EmptySequence reason lets upstream
        /// code distinguish "user gave us nothing" from "the model choked".
        /// </summary>
        [Test]
        public void PredictRetentionTime_EmptyBaseSequence_ReturnsNullWithEmptySequence()
        {
            var model = new TestableRetentionTimeModel();
            var peptide = new StubPeptide { BaseSequence = "", FullSequence = "" };

            var result = model.PredictRetentionTime(peptide, out var reason);

            Assert.That(result, Is.Null);
            Assert.That(reason, Is.EqualTo(RetentionTimeFailureReason.EmptySequence));
        }

        /// <summary>
        /// Happy path: the batch returns a value keyed by FullSequence, and the
        /// single-peptide method unwraps it. Confirms the list-of-one delegation is
        /// wired up correctly end-to-end.
        /// </summary>
        [Test]
        public void PredictRetentionTime_Success_ReturnsPredictedValue()
        {
            var model = new TestableRetentionTimeModel
            {
                PredictStub = inputs => inputs.Select(i => Pred(i.FullSequence, 42.5)).ToList()
            };
            var peptide = new StubPeptide { BaseSequence = "PEPTIDE", FullSequence = "PEPTIDE" };

            var result = model.PredictRetentionTime(peptide, out var reason);

            Assert.That(result, Is.EqualTo(42.5));
            Assert.That(reason, Is.Null);
        }

        /// <summary>
        /// When the batch produces a null value for the requested FullSequence (model
        /// returned no prediction, or the peptide was filtered), the single-peptide
        /// method must surface PredictionError rather than returning null silently.
        /// </summary>
        [Test]
        public void PredictRetentionTime_BatchReturnsNull_ReturnsNullWithPredictionError()
        {
            var model = new TestableRetentionTimeModel
            {
                PredictStub = inputs => inputs.Select(i => Pred(i.FullSequence, null)).ToList()
            };
            var peptide = new StubPeptide { BaseSequence = "PEPTIDE", FullSequence = "PEPTIDE" };

            var result = model.PredictRetentionTime(peptide, out var reason);

            Assert.That(result, Is.Null);
            Assert.That(reason, Is.EqualTo(RetentionTimeFailureReason.PredictionError));
        }

        /// <summary>
        /// A runtime exception inside the batch call (HTTP failure, deserialization
        /// error, etc.) must be swallowed into PredictionError. Propagating exceptions
        /// from a single-peptide call would crash PEP training loops.
        /// </summary>
        [Test]
        public void PredictRetentionTime_BatchThrows_ReturnsNullWithPredictionError()
        {
            var model = new TestableRetentionTimeModel { ThrowOnPredict = true };
            var peptide = new StubPeptide { BaseSequence = "PEPTIDE", FullSequence = "PEPTIDE" };

            var result = model.PredictRetentionTime(peptide, out var reason);

            Assert.That(result, Is.Null);
            Assert.That(reason, Is.EqualTo(RetentionTimeFailureReason.PredictionError));
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // PredictRetentionTimes (batch)
        // ═══════════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Passing null must throw ArgumentNullException (documented contract).
        /// Silently returning an empty dictionary would mask bugs in callers.
        /// </summary>
        [Test]
        public void PredictRetentionTimes_NullInput_ThrowsArgumentNullException()
        {
            var model = new TestableRetentionTimeModel();

            Assert.That(() => model.PredictRetentionTimes(null!),
                Throws.ArgumentNullException);
        }

        /// <summary>
        /// Empty input returns an empty read-only dictionary and never touches the
        /// HTTP layer. Important for the MetaMorpheus code path that filters by
        /// confidence before calling — an empty filter result should be a no-op.
        /// </summary>
        [Test]
        public void PredictRetentionTimes_EmptyInput_ReturnsEmptyReadOnlyDictionary()
        {
            var model = new TestableRetentionTimeModel
            {
                PredictStub = _ => throw new InvalidOperationException("Predict should not be called for empty input")
            };

            var results = model.PredictRetentionTimes(new List<IRetentionPredictable>());

            Assert.That(results, Is.Empty);
            Assert.That(results, Is.InstanceOf<IReadOnlyDictionary<string, double?>>());
        }

        /// <summary>
        /// Input is deduplicated on FullSequence before the API call. Koina bills per
        /// unique sequence, and the base class promises callers that duplicates are
        /// not re-predicted redundantly at the HTTP layer.
        /// </summary>
        [Test]
        public void PredictRetentionTimes_DuplicateFullSequences_DeduplicatedBeforePredict()
        {
            int predictCallCount = 0;
            var model = new TestableRetentionTimeModel
            {
                PredictStub = inputs =>
                {
                    predictCallCount = inputs.Count;
                    return inputs.Select(i => Pred(i.FullSequence, 10.0)).ToList();
                }
            };

            var peptides = new List<IRetentionPredictable>
            {
                new StubPeptide { BaseSequence = "PEPTIDE", FullSequence = "PEPTIDE" },
                new StubPeptide { BaseSequence = "PEPTIDE", FullSequence = "PEPTIDE" },
                new StubPeptide { BaseSequence = "PROTEIN", FullSequence = "PROTEIN" }
            };

            var results = model.PredictRetentionTimes(peptides);

            Assert.That(predictCallCount, Is.EqualTo(2),
                "Duplicate FullSequence entries must be collapsed before hitting the API");
            Assert.That(results.Count, Is.EqualTo(2));
            Assert.That(results["PEPTIDE"], Is.EqualTo(10.0));
            Assert.That(results["PROTEIN"], Is.EqualTo(10.0));
        }

        /// <summary>
        /// Happy path: all predicted values flow through to the returned dictionary
        /// keyed by FullSequence. This is the common case MetaMorpheus sees during
        /// PEP training.
        /// </summary>
        [Test]
        public void PredictRetentionTimes_Success_ReturnsAllPredictions()
        {
            var model = new TestableRetentionTimeModel
            {
                PredictStub = inputs => new List<PeptideRTPrediction>
                {
                    Pred("PEPTIDE", 25.3),
                    Pred("PROTEIN", 42.1)
                }
            };

            var peptides = new List<IRetentionPredictable>
            {
                new StubPeptide { FullSequence = "PEPTIDE" },
                new StubPeptide { FullSequence = "PROTEIN" }
            };

            var results = model.PredictRetentionTimes(peptides);

            Assert.That(results.Count, Is.EqualTo(2));
            Assert.That(results["PEPTIDE"], Is.EqualTo(25.3));
            Assert.That(results["PROTEIN"], Is.EqualTo(42.1));
        }

        /// <summary>
        /// Defensive fill: if the Predict call returns fewer results than requested
        /// (e.g. the underlying model silently drops some), every unmatched input
        /// must still appear in the output dictionary as null. Callers iterate the
        /// dictionary expecting a complete keyset.
        /// </summary>
        [Test]
        public void PredictRetentionTimes_PredictReturnsFewerResults_MissingKeysFilledWithNull()
        {
            var model = new TestableRetentionTimeModel
            {
                PredictStub = _ => new List<PeptideRTPrediction>
                {
                    Pred("PEPTIDE", 25.3)
                    // PROTEIN intentionally missing
                }
            };

            var peptides = new List<IRetentionPredictable>
            {
                new StubPeptide { FullSequence = "PEPTIDE" },
                new StubPeptide { FullSequence = "PROTEIN" }
            };

            var results = model.PredictRetentionTimes(peptides);

            Assert.That(results.Count, Is.EqualTo(2));
            Assert.That(results["PEPTIDE"], Is.EqualTo(25.3));
            Assert.That(results["PROTEIN"], Is.Null);
        }

        /// <summary>
        /// Catastrophic failure in the Predict layer (HTTP unreachable, bad response,
        /// etc.) must not bubble up. Every input gets a null entry so downstream PEP
        /// code sees sentinel values uniformly rather than half a dictionary plus an
        /// exception.
        /// </summary>
        [Test]
        public void PredictRetentionTimes_PredictThrows_ReturnsAllNullEntries()
        {
            var model = new TestableRetentionTimeModel { ThrowOnPredict = true };

            var peptides = new List<IRetentionPredictable>
            {
                new StubPeptide { FullSequence = "PEPTIDE" },
                new StubPeptide { FullSequence = "PROTEIN" }
            };

            var results = model.PredictRetentionTimes(peptides);

            Assert.That(results.Count, Is.EqualTo(2));
            Assert.That(results["PEPTIDE"], Is.Null);
            Assert.That(results["PROTEIN"], Is.Null);
        }

        /// <summary>
        /// The returned dictionary must be truly immutable — callers cannot cast it
        /// back to Dictionary and add entries. This contract allows the base class to
        /// safely return its internal mapping without defensive copying elsewhere.
        /// </summary>
        [Test]
        public void PredictRetentionTimes_ReturnedDictionary_IsImmutable()
        {
            var model = new TestableRetentionTimeModel
            {
                PredictStub = _ => new List<PeptideRTPrediction> { Pred("PEPTIDE", 1.0) }
            };

            var results = model.PredictRetentionTimes(new List<IRetentionPredictable>
            {
                new StubPeptide { FullSequence = "PEPTIDE" }
            });

            Assert.That(results as Dictionary<string, double?>, Is.Null,
                "Returned dictionary must not be castable to a mutable Dictionary");
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // GetFormattedSequence
        // ═══════════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Null peptide must not NRE. This method is used for diagnostics and
        /// caching — an NRE here would destabilize any code that calls it defensively.
        /// The reason is EmptySequence (matching PredictRetentionTime's null-peptide
        /// handling) so callers branching on FailureReason see the same value across
        /// the IRetentionTimePredictor interface methods.
        /// </summary>
        [Test]
        public void GetFormattedSequence_NullPeptide_ReturnsNullWithEmptySequence()
        {
            var model = new TestableRetentionTimeModel();

            var result = model.GetFormattedSequence(null!, out var reason);

            Assert.That(result, Is.Null);
            Assert.That(reason, Is.EqualTo(RetentionTimeFailureReason.EmptySequence));
        }

        /// <summary>
        /// Empty FullSequence is distinguishable from a cleaning failure. The
        /// EmptySequence reason lets callers avoid counting this as a model-side error.
        /// </summary>
        [Test]
        public void GetFormattedSequence_EmptyFullSequence_ReturnsNullWithEmptySequence()
        {
            var model = new TestableRetentionTimeModel();
            var peptide = new StubPeptide { BaseSequence = "PEPTIDE", FullSequence = "" };

            var result = model.GetFormattedSequence(peptide, out var reason);

            Assert.That(result, Is.Null);
            Assert.That(reason, Is.EqualTo(RetentionTimeFailureReason.EmptySequence));
        }

        /// <summary>
        /// Valid peptide: the internal TryCleanSequence path produces a non-null API
        /// sequence. This confirms GetFormattedSequence uses the same conversion
        /// logic as the batch prediction pipeline rather than naively returning
        /// FullSequence verbatim.
        /// </summary>
        [Test]
        public void GetFormattedSequence_ValidPeptide_ReturnsApiSequence()
        {
            var model = new TestableRetentionTimeModel();
            var peptide = new StubPeptide { BaseSequence = "PEPTIDE", FullSequence = "PEPTIDE" };

            var result = model.GetFormattedSequence(peptide, out var reason);

            Assert.That(result, Is.Not.Null);
            Assert.That(result, Is.Not.Empty);
            Assert.That(reason, Is.Null);
        }

        /// <summary>
        /// When TryCleanSequence cannot produce an API sequence (e.g. sequence contains
        /// invalid characters and the handling mode is ReturnNull), GetFormattedSequence
        /// must return null with IncompatibleModifications. This is the failure reason
        /// MetaMorpheus maps to "skip this peptide" during calibration.
        /// </summary>
        [Test]
        public void GetFormattedSequence_UncleanableSequence_ReturnsNullWithIncompatibleModifications()
        {
            var model = new TestableRetentionTimeModel(SequenceConversionHandlingMode.ReturnNull);
            // '*' is not in the allowed amino-acid regex and there's no mod to strip
            var peptide = new StubPeptide { BaseSequence = "PEP*TIDE", FullSequence = "PEP*TIDE" };

            var result = model.GetFormattedSequence(peptide, out var reason);

            Assert.That(result, Is.Null);
            Assert.That(reason, Is.EqualTo(RetentionTimeFailureReason.IncompatibleModifications));
        }

        /// <summary>
        /// ThrowException handling mode causes TryCleanSequence to throw on bad input.
        /// GetFormattedSequence must catch the exception and return null with
        /// PredictionError rather than letting it escape. This protects callers that
        /// happen to use a strict handling mode elsewhere in the pipeline.
        /// </summary>
        [Test]
        public void GetFormattedSequence_TryCleanThrows_ReturnsNullWithPredictionError()
        {
            var model = new TestableRetentionTimeModel(SequenceConversionHandlingMode.ThrowException);
            var peptide = new StubPeptide { BaseSequence = "PEP*TIDE", FullSequence = "PEP*TIDE" };

            var result = model.GetFormattedSequence(peptide, out var reason);

            Assert.That(result, Is.Null);
            Assert.That(reason, Is.EqualTo(RetentionTimeFailureReason.PredictionError));
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // ResponseToPredictions (protected virtual — invoked via test wrapper)
        // ═══════════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Empty requestInputs produces an empty output list and never attempts to
        /// deserialize — responses are ignored entirely in this case.
        /// </summary>
        [Test]
        public void ResponseToPredictions_EmptyRequestInputs_ReturnsEmpty()
        {
            var model = new TestableRetentionTimeModel();

            var result = model.InvokeResponseToPredictions(
                new[] { "ignored" },
                new List<RetentionTimePredictionInput>());

            Assert.That(result, Is.Empty);
        }

        /// <summary>
        /// If any response deserializes to null (literal "null" JSON, corrupted
        /// payload, etc.), the method must throw rather than silently return wrong
        /// predictions. This fails loud — preferred for an HTTP client layer.
        /// </summary>
        [Test]
        public void ResponseToPredictions_NullDeserialization_Throws()
        {
            var model = new TestableRetentionTimeModel();
            var requests = new List<RetentionTimePredictionInput>
            {
                new("PEPTIDE") { ValidatedFullSequence = "PEPTIDE" }
            };

            Assert.That(() => model.InvokeResponseToPredictions(new[] { "null" }, requests),
                Throws.Exception.With.Message.Contains("deserialization"));
        }

        /// <summary>
        /// Response count must match input count. A mismatch indicates the server
        /// silently dropped (or added) entries and the results cannot be safely
        /// mapped back to the requesting peptides.
        /// </summary>
        [Test]
        public void ResponseToPredictions_CountMismatch_Throws()
        {
            var model = new TestableRetentionTimeModel();
            var requests = new List<RetentionTimePredictionInput>
            {
                new("PEPTIDE") { ValidatedFullSequence = "PEPTIDE" },
                new("PROTEIN") { ValidatedFullSequence = "PROTEIN" }
            };
            // Response contains only ONE data point for TWO inputs
            var response = BuildJsonResponse(42.5);

            Assert.That(() => model.InvokeResponseToPredictions(new[] { response }, requests),
                Throws.Exception.With.Message.Contains("number of predictions"));
        }

        /// <summary>
        /// Regression test for the pre-Rev-2 index bug: when inputs had been filtered
        /// by validation, the mapping used the wrong ValidatedFullSequence. The fix
        /// changed `ModelInputs[i]` → `requestInputs[i]`, and this test pins that
        /// contract so the bug cannot silently return.
        ///
        /// We also verify that IsIndexed propagates from the model's
        /// IsIndexedRetentionTimeModel setting, and that the Warning field carries
        /// through from the SequenceWarning on the input.
        /// </summary>
        [Test]
        public void ResponseToPredictions_Success_MapsValidatedFullSequenceFromRequestInputs()
        {
            var model = new TestableRetentionTimeModel();
            var requests = new List<RetentionTimePredictionInput>
            {
                new("PEPTIDE") { ValidatedFullSequence = "PEPTIDE_api", SequenceWarning = null },
                new("PROTEIN") { ValidatedFullSequence = "PROTEIN_api", SequenceWarning = new WarningException("minor") }
            };
            var response = BuildJsonResponse(42.5, 17.25);

            var result = model.InvokeResponseToPredictions(new[] { response }, requests);

            Assert.That(result.Count, Is.EqualTo(2));

            Assert.That(result[0].FullSequence, Is.EqualTo("PEPTIDE"));
            Assert.That(result[0].ValidatedFullSequence, Is.EqualTo("PEPTIDE_api"));
            Assert.That(result[0].PredictedRetentionTime, Is.EqualTo(42.5));
            Assert.That(result[0].IsIndexed, Is.True, "IsIndexed must reflect the model's IsIndexedRetentionTimeModel flag");
            Assert.That(result[0].Warning, Is.Null);

            Assert.That(result[1].FullSequence, Is.EqualTo("PROTEIN"));
            Assert.That(result[1].ValidatedFullSequence, Is.EqualTo("PROTEIN_api"));
            Assert.That(result[1].PredictedRetentionTime, Is.EqualTo(17.25));
            Assert.That(result[1].Warning, Is.Not.Null);
            Assert.That(result[1].Warning!.Message, Is.EqualTo("minor"));
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // AsyncThrottledPredictor / Predict — paths that don't require HTTP
        // ═══════════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Null input to Predict() must not NRE. The base class short-circuits to an
        /// empty predictions list. This is the contract MetaMorpheus relies on when
        /// it calls Predict with a filter result that may have produced zero entries.
        /// </summary>
        [Test]
        public void Predict_NullInputs_ReturnsEmptyAndLeavesPredictionsEmpty()
        {
            var model = new TestableRetentionTimeModel { BypassHttpWithStub = false };

            var result = model.Predict(null!);

            Assert.That(result, Is.Empty);
            Assert.That(model.Predictions, Is.Empty);
        }

        /// <summary>
        /// Empty input list is indistinguishable from null for this contract — both
        /// short-circuit before HTTP.
        /// </summary>
        [Test]
        public void Predict_EmptyInputs_ReturnsEmpty()
        {
            var model = new TestableRetentionTimeModel { BypassHttpWithStub = false };

            var result = model.Predict(new List<RetentionTimePredictionInput>());

            Assert.That(result, Is.Empty);
        }

        /// <summary>
        /// If every input fails validation, no HTTP call is made and the realignment
        /// loop produces a placeholder PeptideRTPrediction for each input with:
        ///   • PredictedRetentionTime = null
        ///   • Warning != null (either the sequence-specific warning from TryCleanSequence
        ///     or the fallback "Input was invalid and skipped during prediction" message)
        /// This is the path that exercises the `else` branch of the realignment loop
        /// without needing to stub HTTP.
        /// </summary>
        [Test]
        public void Predict_AllInputsInvalid_ReturnsPlaceholdersWithWarnings()
        {
            var model = new TestableRetentionTimeModel(SequenceConversionHandlingMode.ReturnNull)
            {
                BypassHttpWithStub = false
            };
            var inputs = new List<RetentionTimePredictionInput>
            {
                new("PEP*TIDE"),    // invalid character
                new("SEL2CYS")      // '2' is not a valid AA
            };

            var result = model.Predict(inputs);

            Assert.That(result.Count, Is.EqualTo(2));
            Assert.That(result.All(r => r.PredictedRetentionTime is null), Is.True);
            Assert.That(result.All(r => r.Warning is not null), Is.True,
                "Every invalid input must carry a warning explaining why it was skipped");
            Assert.That(model.ValidInputsMask, Is.All.False);
        }

        // ── PredictRetentionTimeEquivalents: IRetentionTimePredictor batch contract ──

        [Test]
        public void PredictRetentionTimeEquivalents_NullPeptides_ThrowsArgumentNullException()
        {
            // Pins the input-null guard at the top of PredictRetentionTimeEquivalents
            // so a typo'd or unset peptide list fails fast instead of NRE'ing deeper in.
            var model = new TestableRetentionTimeModel();

            Assert.That(
                () => model.PredictRetentionTimeEquivalents(null!),
                Throws.TypeOf<ArgumentNullException>()
                      .With.Property("ParamName").EqualTo("peptides"));
        }

        [Test]
        public void PredictRetentionTimeEquivalents_SuccessfulBatch_ReturnsTuplesKeyedByOriginalPeptides()
        {
            // Pins the happy-path reshape: each input peptide produces one tuple in the
            // output, the PredictedValue carries through from the batch result, and the
            // FailureReason is null. The original IRetentionPredictable is round-tripped
            // through .Peptide so callers can join results back to their input objects.
            var model = new TestableRetentionTimeModel
            {
                PredictStub = inputs => inputs
                    .Select(inp => Pred(inp.FullSequence, 12.5))
                    .ToList()
            };
            var peptides = new[]
            {
                new StubPeptide { FullSequence = "PEPTIDE",  BaseSequence = "PEPTIDE"  },
                new StubPeptide { FullSequence = "ANALYSIS", BaseSequence = "ANALYSIS" },
            };

            var result = model.PredictRetentionTimeEquivalents(peptides);

            Assert.That(result.Count, Is.EqualTo(2));
            Assert.That(result.All(r => r.PredictedValue.HasValue), Is.True);
            Assert.That(result.All(r => r.PredictedValue == 12.5), Is.True);
            Assert.That(result.All(r => r.FailureReason is null), Is.True);
            Assert.That(result.Select(r => r.Peptide.FullSequence).OrderBy(s => s),
                Is.EqualTo(new[] { "ANALYSIS", "PEPTIDE" }).AsCollection,
                "Each input peptide must round-trip through the .Peptide tuple slot");
        }

        [Test]
        public void PredictRetentionTimeEquivalents_FailedPrediction_PreservesWarningAsIncompatibleMods()
        {
            // Pins the failure-reshape branch: when the batch result has no value for an
            // input peptide AND the matching Predictions entry carries a Warning, the
            // tuple's FailureReason is IncompatibleModifications (not the catch-all
            // PredictionError). Lets callers distinguish "model bug / network error"
            // from "this peptide had a sequence/mod the model couldn't handle."
            var model = new TestableRetentionTimeModel
            {
                PredictStub = inputs => inputs.Select(inp =>
                    new PeptideRTPrediction(
                        FullSequence:           inp.FullSequence,
                        ValidatedFullSequence:  inp.FullSequence,
                        PredictedRetentionTime: null,
                        IsIndexed:              true,
                        Warning:                new WarningException("Mod not in UNIMOD"))
                ).ToList()
            };
            var peptides = new[]
            {
                new StubPeptide { FullSequence = "PEPTIDE", BaseSequence = "PEPTIDE" },
            };

            var result = model.PredictRetentionTimeEquivalents(peptides);

            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].PredictedValue, Is.Null);
            Assert.That(result[0].FailureReason,
                Is.EqualTo(RetentionTimeFailureReason.IncompatibleModifications),
                "A null prediction with a non-null Warning must map to IncompatibleModifications");
        }
    }
}
