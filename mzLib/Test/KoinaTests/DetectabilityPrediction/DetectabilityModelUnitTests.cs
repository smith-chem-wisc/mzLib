using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;

namespace Test.KoinaTests.DetectabilityPrediction
{
    /// <summary>
    /// Unit tests for DetectabilityModel that run without any network access.
    /// Tests ResponseToPredictions logic using a test double.
    /// </summary>
    [TestFixture]
    public class DetectabilityModelUnitTests
    {
        [Test]
        public void ResponseToPredictions_ParsesFourClassProbabilities()
        {
            var model = new TestableDetectabilityModel();
            var inputs = new List<DetectabilityPredictionInput>
            {
                new("PEPTIDEK") { ValidatedFullSequence = "PEPTIDEK" }
            };

            var response = BuildJsonResponse(new[] { 0.1, 0.2, 0.3, 0.4 });
            var predictions = model.TestResponseToPredictions(new[] { response }, inputs);

            Assert.That(predictions.Count, Is.EqualTo(1));
            var probs = predictions[0].DetectabilityProbabilities!.Value;
            Assert.That(probs.NotDetectable, Is.EqualTo(0.1));
            Assert.That(probs.LowDetectability, Is.EqualTo(0.2));
            Assert.That(probs.IntermediateDetectability, Is.EqualTo(0.3));
            Assert.That(probs.HighDetectability, Is.EqualTo(0.4));
        }

        [Test]
        public void ResponseToPredictions_ParsesMultiplePeptides()
        {
            var model = new TestableDetectabilityModel();
            var inputs = new List<DetectabilityPredictionInput>
            {
                new("PEPTIDEK") { ValidatedFullSequence = "PEPTIDEK" },
                new("ANOTHERSEQ") { ValidatedFullSequence = "ANOTHERSEQ" }
            };

            var response = BuildJsonResponse(new[] { 0.1, 0.2, 0.3, 0.4, 0.2, 0.3, 0.1, 0.4 });
            var predictions = model.TestResponseToPredictions(new[] { response }, inputs);

            Assert.That(predictions.Count, Is.EqualTo(2));

            var probs1 = predictions[0].DetectabilityProbabilities!.Value;
            Assert.That(probs1.NotDetectable, Is.EqualTo(0.1));

            var probs2 = predictions[1].DetectabilityProbabilities!.Value;
            Assert.That(probs2.NotDetectable, Is.EqualTo(0.2));
        }

        [Test]
        public void ResponseToPredictions_PreservesSequenceWarnings()
        {
            var model = new TestableDetectabilityModel();
            var inputs = new List<DetectabilityPredictionInput>
            {
                new("PEPTIDEK") { ValidatedFullSequence = "PEPTIDEK", SequenceWarning = new System.ComponentModel.WarningException("mod removed") }
            };

            var response = BuildJsonResponse(new[] { 0.1, 0.2, 0.3, 0.4 });
            var predictions = model.TestResponseToPredictions(new[] { response }, inputs);

            Assert.That(predictions[0].Warning, Is.Not.Null);
            Assert.That(predictions[0].Warning!.Message, Is.EqualTo("mod removed"));
        }

        [Test]
        public void ResponseToPredictions_ReturnsEmptyForEmptyInputs()
        {
            var model = new TestableDetectabilityModel();
            var response = BuildJsonResponse(new[] { 0.1, 0.2, 0.3, 0.4 });
            var predictions = model.TestResponseToPredictions(new[] { response }, new List<DetectabilityPredictionInput>());

            Assert.That(predictions, Is.Empty);
        }

        [Test]
        public void ResponseToPredictions_ThrowsOnDeserializationFailure()
        {
            var model = new TestableDetectabilityModel();
            var inputs = new List<DetectabilityPredictionInput>
            {
                new("PEPTIDEK") { ValidatedFullSequence = "PEPTIDEK" }
            };

            Assert.Throws<Exception>(() => model.TestResponseToPredictions(new[] { "null" }, inputs));
        }

        private static string BuildJsonResponse(double[] probabilities)
        {
            var dataJson = string.Join(",", probabilities.Select(v => v.ToString(System.Globalization.CultureInfo.InvariantCulture)));
            return $"{{\"outputs\":[{{\"name\":\"detectability\",\"datatype\":\"FP32\",\"shape\":[{probabilities.Length}],\"data\":[{dataJson}]}}]}}";
        }

        private sealed class TestableDetectabilityModel : DetectabilityModel
        {
            private static readonly ISequenceConverter Converter = CreateUnimodConverter(
                UnimodSequenceFormatSchema.Instance, new HashSet<int>());

            public override string ModelName => "TestDetectabilityModel";
            public override int MaxBatchSize => 128;
            public override int MaxNumberOfBatchesPerRequest { get; init; } = 1;
            public override int ThrottlingDelayInMilliseconds { get; init; } = 0;
            public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 0;
            public override int MaxPeptideLength => 40;
            public override int MinPeptideLength => 1;
            public override int NumberOfDetectabilityClasses => 4;
            public override List<string> DetectabilityClasses => new() { "Not Detectable", "Low Detectability", "Intermediate Detectability", "High Detectability" };
            public override SequenceConversionHandlingMode ModHandlingMode { get; init; }

            public TestableDetectabilityModel() : base(Converter) { }

            protected override List<Dictionary<string, object>> ToBatchedRequests(List<DetectabilityPredictionInput> validInputs)
                => new();

            public List<PeptideDetectabilityPrediction> TestResponseToPredictions(
                IReadOnlyList<string> responses, List<DetectabilityPredictionInput> requestInputs)
            {
                ModelInputs = requestInputs;
                return ResponseToPredictions(responses, requestInputs);
            }
        }
    }
}
