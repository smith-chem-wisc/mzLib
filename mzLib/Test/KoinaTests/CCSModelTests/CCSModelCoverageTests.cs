using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.CCSModels;
using PredictionClients.Koina.Util;

namespace Test.KoinaTests.CCSModelTests
{
    /// <summary>
    /// No-network coverage for CollisionalCrossSectionModel validation, the missing-output
    /// error path, and the all-invalid Predict realignment.
    /// </summary>
    [TestFixture]
    public class CCSModelCoverageTests
    {
        [Test]
        public void Validate_RejectsNonPhysicalCharge()
        {
            var model = new TestableCcsModel();
            Assert.That(model.TestValidate(new CCSPredictionInput("PEPTIDEK", 0), out var w0), Is.False);
            Assert.That(w0!.Message, Does.Contain("not a valid physical charge"));
            Assert.That(model.TestValidate(new CCSPredictionInput("PEPTIDEK", -2), out _), Is.False);
        }

        [Test]
        public void Validate_RejectsChargeOutsideAllowedSet()
        {
            var model = new TestableCcsModel();
            Assert.That(model.TestValidate(new CCSPredictionInput("PEPTIDEK", 9), out var w), Is.False);
            Assert.That(w!.Message, Does.Contain("Precursor charge"));
        }

        [Test]
        public void Validate_AcceptsAllowedCharge()
        {
            var model = new TestableCcsModel();
            Assert.That(model.TestValidate(new CCSPredictionInput("PEPTIDEK", 2), out var w), Is.True);
            Assert.That(w, Is.Null);
        }

        [Test]
        public void Validate_NullCharges_SkipsSetCheckButKeepsPhysicalGuard()
        {
            var model = new TestableCcsModel(nullCharges: true);
            // null = charge not applicable: a charge outside any set is accepted...
            Assert.That(model.TestValidate(new CCSPredictionInput("PEPTIDEK", 9), out var w9), Is.True);
            Assert.That(w9, Is.Null);
            // ...but the physical charge >= 1 guard still applies.
            Assert.That(model.TestValidate(new CCSPredictionInput("PEPTIDEK", 0), out _), Is.False);
        }

        [Test]
        public void Validate_ThrowMode_ThrowsOnNonPhysicalCharge()
        {
            var model = new TestableCcsModel(IncompatibleParameterHandlingMode.ThrowException);
            Assert.Throws<ArgumentException>(() => model.TestValidate(new CCSPredictionInput("PEPTIDEK", 0), out _));
        }

        [Test]
        public void ResponseToPredictions_ThrowsWhenCcsOutputMissing()
        {
            var model = new TestableCcsModel();
            var inputs = new List<CCSPredictionInput> { new("PEPTIDEK", 2) { ValidatedFullSequence = "PEPTIDEK" } };
            var response = "{\"outputs\":[{\"name\":\"not_ccs\",\"datatype\":\"FP32\",\"shape\":[1],\"data\":[250.0]}]}";

            Assert.That(() => model.TestResponseToPredictions(new[] { response }, inputs),
                Throws.Exception.With.Message.Contains("ccs"));
        }

        [Test]
        public void Predict_AllInputsInvalid_ReturnsPlaceholdersWithWarnings()
        {
            var model = new IM2Deep();
            var inputs = new List<CCSPredictionInput>
            {
                new("PEPTIDEK", 0),     // non-physical charge
                new("PEP*TIDE", 2)      // invalid sequence
            };

            var predictions = model.Predict(inputs);

            Assert.That(predictions.Count, Is.EqualTo(2));
            Assert.That(predictions.All(p => p.PredictedCCS is null), Is.True);
            Assert.That(predictions.All(p => p.Warning is not null), Is.True);
            Assert.That(model.ValidInputsMask, Is.All.False);
        }

        [Test]
        public void Predict_EmptyInput_ReturnsEmpty()
        {
            Assert.That(new IM2Deep().Predict(new List<CCSPredictionInput>()), Is.Empty);
        }

        private sealed class TestableCcsModel : CollisionalCrossSectionModel
        {
            private static readonly ISequenceConverter Conv = CreateUnimodConverter(
                UnimodSequenceFormatSchema.Instance, new HashSet<int> { 35, 4 });

            public override string ModelName => "TestCCSModel";
            public override int MaxBatchSize => 1000;
            public override int MaxNumberOfBatchesPerRequest { get; init; } = 1;
            public override int ThrottlingDelayInMilliseconds { get; init; } = 0;
            public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 0;
            public override int MaxPeptideLength => 50;
            public override int MinPeptideLength => 1;
            public override SequenceConversionHandlingMode ModHandlingMode { get; init; }
            public override HashSet<int>? AllowedPrecursorCharges { get; }

            public TestableCcsModel(
                IncompatibleParameterHandlingMode paramMode = IncompatibleParameterHandlingMode.ReturnNull,
                bool nullCharges = false)
                : base(Conv)
            {
                ParameterHandlingMode = paramMode;
                AllowedPrecursorCharges = nullCharges ? null : new HashSet<int> { 1, 2, 3, 4, 5, 6 };
            }

            protected override List<Dictionary<string, object>> ToBatchedRequests(List<CCSPredictionInput> validInputs)
                => new();

            public bool TestValidate(CCSPredictionInput input, out System.ComponentModel.WarningException? warning)
                => ValidateModelSpecificInputs(input, out warning);

            public List<PeptideCCSPrediction> TestResponseToPredictions(
                IReadOnlyList<string> responses, List<CCSPredictionInput> requestInputs)
                => ResponseToPredictions(responses, requestInputs);
        }
    }
}
