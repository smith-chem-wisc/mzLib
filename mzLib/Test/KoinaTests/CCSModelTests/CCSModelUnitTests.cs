using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;

namespace Test.KoinaTests.CCSModelTests
{
    /// <summary>
    /// Unit tests for CollisionalCrossSectionModel that run without any network access.
    /// Tests ResponseToPredictions logic using a test double.
    /// </summary>
    [TestFixture]
    public class CCSModelUnitTests
    {
        [Test]
        public void ResponseToPredictions_ParsesSingleCCSValuePerPeptide()
        {
            var model = new TestableCCSModel();
            var inputs = new List<CCSPredictionInput>
            {
                new("PEPTIDEK", 2) { ValidatedFullSequence = "PEPTIDEK" },
                new("ANOTHERSEQ", 3) { ValidatedFullSequence = "ANOTHERSEQ" }
            };

            var response = BuildJsonResponse(new[] { 250.5, 310.2 });
            var predictions = model.TestResponseToPredictions(new[] { response }, inputs);

            Assert.That(predictions.Count, Is.EqualTo(2));
            Assert.That(predictions[0].PredictedCCS, Is.EqualTo(250.5));
            Assert.That(predictions[1].PredictedCCS, Is.EqualTo(310.2));
        }

        [Test]
        public void ResponseToPredictions_ThrowsWhenCountMismatch()
        {
            var model = new TestableCCSModel();
            var inputs = new List<CCSPredictionInput>
            {
                new("PEPTIDEK", 2) { ValidatedFullSequence = "PEPTIDEK" },
                new("ANOTHERSEQ", 3) { ValidatedFullSequence = "ANOTHERSEQ" }
            };

            var response = BuildJsonResponse(new[] { 250.5 });
            Assert.Throws<Exception>(() => model.TestResponseToPredictions(new[] { response }, inputs));
        }

        [Test]
        public void ResponseToPredictions_PreservesSequenceWarnings()
        {
            var model = new TestableCCSModel();
            var inputs = new List<CCSPredictionInput>
            {
                new("PEPTIDEK", 2) { ValidatedFullSequence = "PEPTIDEK", SequenceWarning = new System.ComponentModel.WarningException("minor warning") }
            };

            var response = BuildJsonResponse(new[] { 250.5 });
            var predictions = model.TestResponseToPredictions(new[] { response }, inputs);

            Assert.That(predictions[0].Warning, Is.Not.Null);
            Assert.That(predictions[0].Warning!.Message, Is.EqualTo("minor warning"));
        }

        [Test]
        public void ResponseToPredictions_ReturnsEmptyForEmptyInputs()
        {
            var model = new TestableCCSModel();
            var response = BuildJsonResponse(new[] { 250.5 });
            var predictions = model.TestResponseToPredictions(new[] { response }, new List<CCSPredictionInput>());

            Assert.That(predictions, Is.Empty);
        }

        private static string BuildJsonResponse(double[] ccsValues)
        {
            var dataJson = string.Join(",", ccsValues.Select(v => v.ToString(System.Globalization.CultureInfo.InvariantCulture)));
            return $"{{\"outputs\":[{{\"name\":\"ccs\",\"datatype\":\"FP32\",\"shape\":[{ccsValues.Length}],\"data\":[{dataJson}]}}]}}";
        }

        private sealed class TestableCCSModel : CollisionalCrossSectionModel
        {
            private static readonly IReadOnlySet<int> SupportedUnimodIds = new HashSet<int> { 35, 4 };
            private static readonly ISequenceConverter Converter = CreateUnimodConverter(
                UnimodSequenceFormatSchema.Instance, SupportedUnimodIds);

            public override string ModelName => "TestCCSModel";
            public override int MaxBatchSize => 1000;
            public override int MaxNumberOfBatchesPerRequest { get; init; } = 1;
            public override int ThrottlingDelayInMilliseconds { get; init; } = 0;
            public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 0;
            public override int MaxPeptideLength => 50;
            public override int MinPeptideLength => 1;
            public override HashSet<int>? AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6 };
            public override SequenceConversionHandlingMode ModHandlingMode { get; init; }

            public TestableCCSModel() : base(Converter) { }

            protected override List<Dictionary<string, object>> ToBatchedRequests(List<CCSPredictionInput> validInputs)
                => new();

            public List<PeptideCCSPrediction> TestResponseToPredictions(
                IReadOnlyList<string> responses, List<CCSPredictionInput> requestInputs)
                => ResponseToPredictions(responses, requestInputs);
        }
    }
}
