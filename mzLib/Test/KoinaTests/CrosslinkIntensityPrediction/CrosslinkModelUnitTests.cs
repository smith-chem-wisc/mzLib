using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Client;
using PredictionClients.Koina.SupportedModels.CrosslinkIntensityModels;
using PredictionClients.Koina.Util;

namespace Test.KoinaTests.CrosslinkIntensityPrediction
{
    /// <summary>
    /// Unit tests for CrosslinkFragmentIntensityModel that run without any network access.
    /// Covers TryCleanSequence (UNIMOD-marker validation), ValidateModelSpecificInputs,
    /// ExtractOutputs, ResponseToPredictions, and the all-invalid realignment path.
    /// </summary>
    [TestFixture]
    public class CrosslinkModelUnitTests
    {
        private static string BuildJsonResponse(string[] annotations, double[] mz, double[] intensities)
        {
            var ann = string.Join(",", annotations.Select(a => $"\"{a}\""));
            var m = string.Join(",", mz.Select(v => v.ToString(System.Globalization.CultureInfo.InvariantCulture)));
            var ints = string.Join(",", intensities.Select(v => v.ToString(System.Globalization.CultureInfo.InvariantCulture)));
            return $"{{\"outputs\":[" +
                   $"{{\"name\":\"annotation\",\"datatype\":\"BYTES\",\"shape\":[{annotations.Length}],\"data\":[{ann}]}}," +
                   $"{{\"name\":\"mz\",\"datatype\":\"FP32\",\"shape\":[{mz.Length}],\"data\":[{m}]}}," +
                   $"{{\"name\":\"intensities\",\"datatype\":\"FP32\",\"shape\":[{intensities.Length}],\"data\":[{ints}]}}]}}";
        }

        // ── TryCleanSequence override ──────────────────────────────────────────────

        [Test]
        public void TryCleanSequence_ValidUnimodMarker_PassesThroughUnchanged()
        {
            var model = new TestablePairedModel();
            var result = model.TestTryCleanSequence("PEPTIDEK[UNIMOD:1896]", out var api, out var warning);

            Assert.That(result, Is.EqualTo("PEPTIDEK[UNIMOD:1896]"));
            Assert.That(api, Is.EqualTo("PEPTIDEK[UNIMOD:1896]"));
            Assert.That(warning, Is.Null);
        }

        [Test]
        public void TryCleanSequence_DisallowedUnimodId_ReturnsNullWithWarning()
        {
            var model = new TestablePairedModel();
            var result = model.TestTryCleanSequence("PEPS[UNIMOD:21]IDEK[UNIMOD:1896]", out var api, out var warning);

            Assert.That(result, Is.Null);
            Assert.That(api, Is.Null);
            Assert.That(warning, Is.Not.Null);
            Assert.That(warning!.Message, Does.Contain("UNIMOD:21"));
        }

        [Test]
        public void TryCleanSequence_NonUnimodNotation_ReturnsNullWithWarning()
        {
            var model = new TestablePairedModel();
            var result = model.TestTryCleanSequence("K[Common Variable:Acetyl on K]PEPTIDEK", out var api, out var warning);

            Assert.That(result, Is.Null);
            Assert.That(warning, Is.Not.Null);
        }

        [Test]
        public void TryCleanSequence_InvalidBaseSequence_ReturnsNullWithWarning()
        {
            var model = new TestablePairedModel();
            var result = model.TestTryCleanSequence("PEP*TIDE", out var api, out var warning);

            Assert.That(result, Is.Null);
            Assert.That(warning, Is.Not.Null);
        }

        [Test]
        public void TryCleanSequence_ThrowMode_ThrowsOnInvalidBaseSequence()
        {
            var model = new TestablePairedModel(SequenceConversionHandlingMode.ThrowException);
            Assert.Throws<ArgumentException>(() => model.TestTryCleanSequence("PEP*TIDE", out _, out _));
        }

        [Test]
        public void TryCleanSequence_ThrowMode_ThrowsOnDisallowedUnimodId()
        {
            var model = new TestablePairedModel(SequenceConversionHandlingMode.ThrowException);
            Assert.Throws<ArgumentException>(
                () => model.TestTryCleanSequence("PEPS[UNIMOD:21]IDEK[UNIMOD:1896]", out _, out _));
        }

        [Test]
        public void RequiresBetaSequence_TrueForPairedModels_FalseForSingleSequence()
        {
            Assert.That(new Prosit2023IntensityXLCMS3().RequiresBetaSequence, Is.False);
            Assert.That(new Prosit2023IntensityXLCMS2().RequiresBetaSequence, Is.True);
            Assert.That(new Prosit2024IntensityXLNMS2().RequiresBetaSequence, Is.True);
        }

        // ── ValidateModelSpecificInputs ────────────────────────────────────────────

        [Test]
        public void Validate_PairedModelRejectsNullBeta()
        {
            var model = new TestablePairedModel();
            var input = new CrosslinkIntensityPredictionInput("PEPTIDEK[UNIMOD:1896]", null, 2, 35);

            Assert.That(model.TestValidate(input, out var warning), Is.False);
            Assert.That(warning!.Message, Does.Contain("Beta sequence is required"));
        }

        [Test]
        public void Validate_RejectsUnsupportedCharge()
        {
            var model = new TestablePairedModel();
            var input = new CrosslinkIntensityPredictionInput("PEPTIDEK[UNIMOD:1896]", "ACDEFGK[UNIMOD:1896]", 9, 35);

            Assert.That(model.TestValidate(input, out var warning), Is.False);
            Assert.That(warning!.Message, Does.Contain("Precursor charge"));
        }

        [Test]
        public void Validate_RejectsMissingCollisionEnergyWhenRequired()
        {
            // Empty AllowedCollisionEnergies means CE is required but any value is accepted.
            var model = new TestablePairedModel();
            var input = new CrosslinkIntensityPredictionInput("PEPTIDEK[UNIMOD:1896]", "ACDEFGK[UNIMOD:1896]", 2, null);

            Assert.That(model.TestValidate(input, out var warning), Is.False);
            Assert.That(warning!.Message, Does.Contain("CollisionEnergy"));
        }

        [Test]
        public void Validate_ThrowMode_ThrowsOnMissingBeta()
        {
            var model = new TestablePairedModel(SequenceConversionHandlingMode.ReturnNull, IncompatibleParameterHandlingMode.ThrowException);
            var input = new CrosslinkIntensityPredictionInput("PEPTIDEK[UNIMOD:1896]", null, 2, 35);

            Assert.Throws<ArgumentException>(() => model.TestValidate(input, out _));
        }

        [Test]
        public void Validate_AllParametersValid_ReturnsTrue()
        {
            var model = new TestablePairedModel();
            var input = new CrosslinkIntensityPredictionInput("PEPTIDEK[UNIMOD:1896]", "ACDEFGK[UNIMOD:1896]", 3, 30);

            Assert.That(model.TestValidate(input, out var warning), Is.True);
            Assert.That(warning, Is.Null);
        }

        // ── ExtractOutputs ─────────────────────────────────────────────────────────

        [Test]
        public void ExtractOutputs_ThrowsWhenMissingOutput()
        {
            var model = new TestablePairedModel();
            var response = new ResponseJSONStruct
            {
                Outputs = new List<OutputJSONStruct>
                {
                    new() { Name = "annotation", Data = new List<object> { "b1+1" } },
                    new() { Name = "mz", Data = new List<object> { 100.0 } }
                }
            };

            Assert.Throws<Exception>(() => model.TestExtractOutputs(response));
        }

        // ── ResponseToPredictions ──────────────────────────────────────────────────

        [Test]
        public void ResponseToPredictions_ParsesAndSkipsImpossibleIons()
        {
            var model = new TestablePairedModel();
            var inputs = new List<CrosslinkIntensityPredictionInput>
            {
                new("PEPTIDEK[UNIMOD:1896]", "ACDEFGK[UNIMOD:1896]", 2, 35)
                {
                    ValidatedAlphaSequence = "PEPTIDEK[UNIMOD:1896]",
                    ValidatedBetaSequence = "ACDEFGK[UNIMOD:1896]"
                }
            };
            var response = BuildJsonResponse(
                annotations: new[] { "b1+1", "y1+1", "a1+1" },
                mz: new[] { 150.0, 200.0, 100.0 },
                intensities: new[] { 0.8, -1.0, 0.6 });

            var predictions = model.TestResponseToPredictions(new[] { response }, inputs);

            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].FragmentIntensities, Is.EquivalentTo(new[] { 0.8, 0.6 }));
            Assert.That(predictions[0].ValidatedAlphaSequence, Is.EqualTo("PEPTIDEK[UNIMOD:1896]"));
            Assert.That(predictions[0].ValidatedBetaSequence, Is.EqualTo("ACDEFGK[UNIMOD:1896]"));
        }

        [Test]
        public void ResponseToPredictions_SplitsFragmentsAcrossMultiplePeptides()
        {
            var model = new TestablePairedModel();
            var inputs = new List<CrosslinkIntensityPredictionInput>
            {
                new("PEPTIDEK[UNIMOD:1896]", "ACDEFGK[UNIMOD:1896]", 2, 35) { ValidatedAlphaSequence = "PEPTIDEK[UNIMOD:1896]", ValidatedBetaSequence = "ACDEFGK[UNIMOD:1896]" },
                new("LMNPEK[UNIMOD:1896]", "PEPK[UNIMOD:1896]", 3, 35) { ValidatedAlphaSequence = "LMNPEK[UNIMOD:1896]", ValidatedBetaSequence = "PEPK[UNIMOD:1896]" }
            };
            var response = BuildJsonResponse(
                annotations: new[] { "b1+1", "y1+1" },
                mz: new[] { 150.0, 200.0 },
                intensities: new[] { 0.8, 0.6 });

            var predictions = model.TestResponseToPredictions(new[] { response }, inputs);

            Assert.That(predictions.Count, Is.EqualTo(2));
            Assert.That(predictions[0].FragmentIntensities, Is.EquivalentTo(new[] { 0.8 }));
            Assert.That(predictions[1].FragmentIntensities, Is.EquivalentTo(new[] { 0.6 }));
        }

        [Test]
        public void ResponseToPredictions_EmptyInputs_ReturnsEmpty()
        {
            var model = new TestablePairedModel();
            var response = BuildJsonResponse(new[] { "b1+1" }, new[] { 150.0 }, new[] { 0.8 });

            Assert.That(model.TestResponseToPredictions(new[] { response }, new List<CrosslinkIntensityPredictionInput>()), Is.Empty);
        }

        // ── Predict realignment (all inputs invalid → no network) ──────────────────

        [Test]
        public void Predict_AllInputsInvalid_ReturnsPlaceholdersWithWarnings()
        {
            var model = new Prosit2023IntensityXLCMS2();
            var inputs = new List<CrosslinkIntensityPredictionInput>
            {
                new("PEPTIDEK[UNIMOD:1896]", null, 2, 35),               // missing beta
                new("PEPS[UNIMOD:21]IDEK[UNIMOD:1896]", "ACDEK[UNIMOD:1896]", 2, 35) // disallowed mod
            };

            var predictions = model.Predict(inputs);

            Assert.That(predictions.Count, Is.EqualTo(2));
            Assert.That(predictions.All(p => p.FragmentAnnotations is null), Is.True);
            Assert.That(predictions.All(p => p.Warning is not null), Is.True);
            Assert.That(model.ValidInputsMask, Is.All.False);
        }

        [Test]
        public void Predict_EmptyInput_ReturnsEmpty()
        {
            var model = new Prosit2023IntensityXLCMS2();
            Assert.That(model.Predict(new List<CrosslinkIntensityPredictionInput>()), Is.Empty);
        }

        // ── Prediction record carries beta sequence info ───────────────────────────

        [Test]
        public void Prediction_Record_PreservesBetaSequence()
        {
            var withBeta = new CrosslinkFragmentIntensityPrediction(
                "ALPHA", "BETA", "ALPHA_v", "BETA_v", 2, null, null, null);
            var singleSeq = new CrosslinkFragmentIntensityPrediction(
                "ALPHA", null, "ALPHA_v", null, 2, null, null, null);

            Assert.That(withBeta.BetaSequence, Is.EqualTo("BETA"));
            Assert.That(singleSeq.BetaSequence, Is.Null);
            Assert.That(singleSeq.ValidatedBetaSequence, Is.Null);
        }

        // ── Test double exposing protected members ─────────────────────────────────

        private sealed class TestablePairedModel : CrosslinkFragmentIntensityModel
        {
            private static readonly IReadOnlySet<int> Ids = new HashSet<int> { 4, 35, 1896 };
            private static readonly ISequenceConverter Conv = CreateUnimodConverter(CrosslinkSchema, Ids);

            public override string ModelName => "TestXLModel";
            public override int MaxBatchSize => 1000;
            public override int MaxNumberOfBatchesPerRequest { get; init; } = 1;
            public override int ThrottlingDelayInMilliseconds { get; init; } = 0;
            public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 0;
            public override int MaxPeptideLength => 30;
            public override int MinPeptideLength => 1;
            public override HashSet<int>? AllowedCollisionEnergies => new HashSet<int>();
            public override IReadOnlySet<int> AllowedUnimodIds => Ids;
            public override SequenceConversionHandlingMode ModHandlingMode { get; init; }
            public override IncompatibleParameterHandlingMode ParameterHandlingMode { get; init; }
            public override int NumberOfPredictedFragmentIons => 348;

            public TestablePairedModel(
                SequenceConversionHandlingMode modMode = SequenceConversionHandlingMode.ReturnNull,
                IncompatibleParameterHandlingMode paramMode = IncompatibleParameterHandlingMode.ReturnNull)
                : base(Conv)
            {
                ModHandlingMode = modMode;
                ParameterHandlingMode = paramMode;
            }

            protected override List<Dictionary<string, object>> ToBatchedRequests(List<CrosslinkIntensityPredictionInput> validInputs)
                => new();

            public string? TestTryCleanSequence(string sequence, out string? api, out System.ComponentModel.WarningException? warning)
                => TryCleanSequence(sequence, out api, out warning);

            public bool TestValidate(CrosslinkIntensityPredictionInput input, out System.ComponentModel.WarningException? warning)
                => ValidateModelSpecificInputs(input, out warning);

            public (List<object>, List<object>, List<object>) TestExtractOutputs(ResponseJSONStruct response)
                => ExtractOutputs(response);

            public List<CrosslinkFragmentIntensityPrediction> TestResponseToPredictions(
                IReadOnlyList<string> responses, List<CrosslinkIntensityPredictionInput> requestInputs)
            {
                ModelInputs = requestInputs;
                return ResponseToPredictions(responses, requestInputs);
            }
        }
    }
}
