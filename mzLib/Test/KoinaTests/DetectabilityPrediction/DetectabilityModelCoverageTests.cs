using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.FlyabilityModels;

namespace Test.KoinaTests.DetectabilityPrediction
{
    /// <summary>
    /// No-network coverage for DetectabilityModel: the all-invalid Predict realignment and a
    /// regression for the validated-sequence index mapping when inputs are filtered.
    /// </summary>
    [TestFixture]
    public class DetectabilityModelCoverageTests
    {
        private static string BuildJsonResponse(double[] probabilities)
        {
            var data = string.Join(",", probabilities.Select(v => v.ToString(System.Globalization.CultureInfo.InvariantCulture)));
            return $"{{\"outputs\":[{{\"name\":\"detectability\",\"datatype\":\"FP32\",\"shape\":[{probabilities.Length}],\"data\":[{data}]}}]}}";
        }

        [Test]
        public void Predict_AllInputsInvalid_ReturnsPlaceholdersWithWarnings()
        {
            var model = new PFly2024FineTuned();
            var inputs = new List<DetectabilityPredictionInput>
            {
                new("PEP*TIDE"),                                   // invalid character
                new(new string('A', 41))                          // exceeds MaxPeptideLength (40)
            };

            var predictions = model.Predict(inputs);

            Assert.That(predictions.Count, Is.EqualTo(2));
            Assert.That(predictions.All(p => p.DetectabilityProbabilities is null), Is.True);
            Assert.That(predictions.All(p => p.Warning is not null), Is.True);
            Assert.That(model.ValidInputsMask, Is.All.False);
        }

        [Test]
        public void Predict_EmptyInput_ReturnsEmpty()
        {
            Assert.That(new PFly2024FineTuned().Predict(new List<DetectabilityPredictionInput>()), Is.Empty);
        }

        [Test]
        public void ResponseToPredictions_MapsValidatedSequenceFromRequestInputs_NotFullModelInputs()
        {
            // When validation filters earlier inputs, requestInputs is a tail subset of
            // ModelInputs. The prediction's ValidatedFullSequence must come from the
            // request input at the same index, not the full ModelInputs list.
            var model = new TestableDetectabilityModel();
            var skipped = new DetectabilityPredictionInput("SKIPPED") { ValidatedFullSequence = null };
            var valid = new DetectabilityPredictionInput("PEPTIDEK") { ValidatedFullSequence = "PEPTIDEK_api" };
            model.SetModelInputs(new List<DetectabilityPredictionInput> { skipped, valid });

            var requestInputs = new List<DetectabilityPredictionInput> { valid };
            var predictions = model.TestResponseToPredictions(new[] { BuildJsonResponse(new[] { 0.1, 0.2, 0.3, 0.4 }) }, requestInputs);

            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].FullSequence, Is.EqualTo("PEPTIDEK"));
            Assert.That(predictions[0].ValidatedFullSequence, Is.EqualTo("PEPTIDEK_api"));
        }

        private sealed class TestableDetectabilityModel : DetectabilityModel
        {
            private static readonly ISequenceConverter Conv = CreateUnimodConverter(
                UnimodSequenceFormatSchema.Instance, new HashSet<int>());

            public override string ModelName => "TestDetectability";
            public override int MaxBatchSize => 128;
            public override int MaxNumberOfBatchesPerRequest { get; init; } = 1;
            public override int ThrottlingDelayInMilliseconds { get; init; } = 0;
            public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 0;
            public override int MaxPeptideLength => 40;
            public override int MinPeptideLength => 1;
            public override int NumberOfDetectabilityClasses => 4;
            public override List<string> DetectabilityClasses => new() { "Not Detectable", "Low Detectability", "Intermediate Detectability", "High Detectability" };
            public override SequenceConversionHandlingMode ModHandlingMode { get; init; }

            public TestableDetectabilityModel() : base(Conv) { }

            protected override List<Dictionary<string, object>> ToBatchedRequests(List<DetectabilityPredictionInput> validInputs)
                => new();

            public void SetModelInputs(List<DetectabilityPredictionInput> modelInputs) => ModelInputs = modelInputs;

            public List<PeptideDetectabilityPrediction> TestResponseToPredictions(
                IReadOnlyList<string> responses, List<DetectabilityPredictionInput> requestInputs)
                => ResponseToPredictions(responses, requestInputs);
        }
    }
}
