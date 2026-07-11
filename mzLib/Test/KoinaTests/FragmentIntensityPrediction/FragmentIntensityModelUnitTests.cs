using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using NUnit.Framework;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Client;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.Util;

namespace Test.KoinaTests.FragmentIntensityPrediction
{
    /// <summary>
    /// Unit tests for FragmentIntensityModel that run without any network access.
    /// Tests ExtractOutputs, ParseFragmentAnnotation, and ResponseToPredictions logic.
    /// </summary>
    [TestFixture]
    public class FragmentIntensityModelUnitTests
    {
        // ═══════════════════════════════════════════════════════════════════════════
        // ExtractOutputs Tests
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        public void ExtractOutputs_FindsOutputsByName_RegardlessOfOrder()
        {
            var response = new ResponseJSONStruct
            {
                Outputs = new List<OutputJSONStruct>
                {
                    new() { Name = "intensities", Data = new List<object> { 0.5, 0.3 } },
                    new() { Name = "mz", Data = new List<object> { 100.0, 200.0 } },
                    new() { Name = "annotation", Data = new List<object> { "b1+1", "y1+1" } }
                }
            };

            var model = new TestableFragmentIntensityModel();
            var (annotations, mz, intensities) = model.TestExtractOutputs(response);

            Assert.That(annotations, Is.EquivalentTo(new[] { "b1+1", "y1+1" }));
            Assert.That(mz, Is.EquivalentTo(new[] { 100.0, 200.0 }));
            Assert.That(intensities, Is.EquivalentTo(new[] { 0.5, 0.3 }));
        }

        [Test]
        public void ExtractOutputs_AcceptsAnnotationsPlural()
        {
            var response = new ResponseJSONStruct
            {
                Outputs = new List<OutputJSONStruct>
                {
                    new() { Name = "annotations", Data = new List<object> { "b1+1" } },
                    new() { Name = "mz", Data = new List<object> { 100.0 } },
                    new() { Name = "intensities", Data = new List<object> { 0.5 } }
                }
            };

            var model = new TestableFragmentIntensityModel();
            var (annotations, mz, intensities) = model.TestExtractOutputs(response);

            Assert.That(annotations, Is.EquivalentTo(new[] { "b1+1" }));
        }

        [Test]
        public void ExtractOutputs_ThrowsWhenMissingOutput()
        {
            var response = new ResponseJSONStruct
            {
                Outputs = new List<OutputJSONStruct>
                {
                    new() { Name = "annotation", Data = new List<object> { "b1+1" } },
                    new() { Name = "mz", Data = new List<object> { 100.0 } }
                }
            };

            var model = new TestableFragmentIntensityModel();
            Assert.Throws<Exception>(() => model.TestExtractOutputs(response));
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // ParseFragmentAnnotation Tests (Default Implementation)
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        public void ParseFragmentAnnotation_Default_ParsesTypePlusCharge()
        {
            var model = new TestableFragmentIntensityModel();

            var result = model.TestParseFragmentAnnotation("b5+1");
            Assert.That(result.FragmentIdentifier, Is.EqualTo("b5"));
            Assert.That(result.Charge, Is.EqualTo(1));

            result = model.TestParseFragmentAnnotation("y10+2");
            Assert.That(result.FragmentIdentifier, Is.EqualTo("y10"));
            Assert.That(result.Charge, Is.EqualTo(2));
        }

        [Test]
        public void ParseFragmentAnnotation_Default_ThrowsOnMissingPlus()
        {
            var model = new TestableFragmentIntensityModel();
            Assert.Throws<Exception>(() => model.TestParseFragmentAnnotation("b5"));
        }

        [Test]
        public void ParseFragmentAnnotation_Default_ThrowsOnInvalidCharge()
        {
            var model = new TestableFragmentIntensityModel();
            Assert.Throws<Exception>(() => model.TestParseFragmentAnnotation("b5+abc"));
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // ParseFragmentAnnotation Tests (Altimeter Override)
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        public void ParseFragmentAnnotation_Altimeter_ParsesCaretCharge()
        {
            var model = new TestableAltimeterModel();

            var result = model.TestParseFragmentAnnotation("y8^2");
            Assert.That(result.FragmentIdentifier, Is.EqualTo("y8"));
            Assert.That(result.Charge, Is.EqualTo(2));

            result = model.TestParseFragmentAnnotation("b3^1");
            Assert.That(result.FragmentIdentifier, Is.EqualTo("b3"));
            Assert.That(result.Charge, Is.EqualTo(1));
        }

        [Test]
        public void ParseFragmentAnnotation_Altimeter_ParsesNoChargeAsCharge1()
        {
            var model = new TestableAltimeterModel();
            var result = model.TestParseFragmentAnnotation("y5");
            Assert.That(result.FragmentIdentifier, Is.EqualTo("y5"));
            Assert.That(result.Charge, Is.EqualTo(1));
        }

        [Test]
        public void ParseFragmentAnnotation_Altimeter_ThrowsOnEmptyAnnotation()
        {
            var model = new TestableAltimeterModel();
            Assert.Throws<Exception>(() => model.TestParseFragmentAnnotation(""));
        }

        [Test]
        public void ParseFragmentAnnotation_Altimeter_ThrowsOnNeutralLossNotation()
        {
            var model = new TestableAltimeterModel();
            Assert.Throws<Exception>(() => model.TestParseFragmentAnnotation("IN+CO"));
        }

        [Test]
        public void ParseFragmentAnnotation_Altimeter_ThrowsOnInvalidCharge()
        {
            var model = new TestableAltimeterModel();
            Assert.Throws<Exception>(() => model.TestParseFragmentAnnotation("y8^abc"));
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // ParseFragmentAnnotation Tests (UniSpec Override)
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        public void ParseFragmentAnnotation_UniSpec_ParsesTypePlusCharge()
        {
            var model = new TestableUniSpecModel();
            var result = model.TestParseFragmentAnnotation("y8+1");
            Assert.That(result.FragmentIdentifier, Is.EqualTo("y8"));
            Assert.That(result.Charge, Is.EqualTo(1));
        }

        [Test]
        public void ParseFragmentAnnotation_UniSpec_ParsesStandardTypeNoCharge()
        {
            var model = new TestableUniSpecModel();
            var result = model.TestParseFragmentAnnotation("a2");
            Assert.That(result.FragmentIdentifier, Is.EqualTo("a2"));
            Assert.That(result.Charge, Is.EqualTo(1));
        }

        [Test]
        public void ParseFragmentAnnotation_UniSpec_ThrowsOnEmptyAnnotation()
        {
            var model = new TestableUniSpecModel();
            Assert.Throws<Exception>(() => model.TestParseFragmentAnnotation(""));
        }

        [Test]
        public void ParseFragmentAnnotation_UniSpec_ThrowsOnInternalFragment()
        {
            var model = new TestableUniSpecModel();
            Assert.Throws<Exception>(() => model.TestParseFragmentAnnotation("Int/EQ/4"));
        }

        [Test]
        public void ParseFragmentAnnotation_UniSpec_ThrowsOnSingleLetterAnnotation()
        {
            var model = new TestableUniSpecModel();
            Assert.Throws<Exception>(() => model.TestParseFragmentAnnotation("IEA"));
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // ParseFragmentAnnotation Tests — Neutral Loss Extraction
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        public void ParseFragmentAnnotation_Default_ExtractsWaterLoss()
        {
            var model = new TestableFragmentIntensityModel();
            var result = model.TestParseFragmentAnnotation("y3-H2O+1");

            Assert.That(result.FragmentIdentifier, Is.EqualTo("y3-H2O"));
            Assert.That(result.Charge, Is.EqualTo(1));
            Assert.That(result.NeutralLossFormula, Is.EqualTo("H2O"));
            Assert.That(result.NeutralLossMass, Is.EqualTo(ChemicalFormula.ParseFormula("H2O").MonoisotopicMass).Within(0.001));
            Assert.That(result.BaseFragmentIdentifier, Is.EqualTo("y3"));
        }

        [Test]
        public void ParseFragmentAnnotation_Default_ExtractsAmmoniaLoss()
        {
            var model = new TestableFragmentIntensityModel();
            var result = model.TestParseFragmentAnnotation("b5-NH3+2");

            Assert.That(result.FragmentIdentifier, Is.EqualTo("b5-NH3"));
            Assert.That(result.Charge, Is.EqualTo(2));
            Assert.That(result.NeutralLossFormula, Is.EqualTo("NH3"));
            Assert.That(result.NeutralLossMass, Is.EqualTo(ChemicalFormula.ParseFormula("NH3").MonoisotopicMass).Within(0.001));
            Assert.That(result.BaseFragmentIdentifier, Is.EqualTo("b5"));
        }

        [Test]
        public void ParseFragmentAnnotation_Default_NoNeutralLoss()
        {
            var model = new TestableFragmentIntensityModel();
            var result = model.TestParseFragmentAnnotation("b5+1");

            Assert.That(result.FragmentIdentifier, Is.EqualTo("b5"));
            Assert.That(result.Charge, Is.EqualTo(1));
            Assert.That(result.NeutralLossFormula, Is.Null);
            Assert.That(result.NeutralLossMass, Is.Null);
            Assert.That(result.BaseFragmentIdentifier, Is.EqualTo("b5"));
        }

        [Test]
        public void ParseFragmentAnnotation_Default_FallsBackOnUnparsableFormula()
        {
            var model = new TestableFragmentIntensityModel();
            var result = model.TestParseFragmentAnnotation("b5-123+1");

            Assert.That(result.FragmentIdentifier, Is.EqualTo("b5-123"));
            Assert.That(result.Charge, Is.EqualTo(1));
            Assert.That(result.NeutralLossFormula, Is.Null);
            Assert.That(result.NeutralLossMass, Is.Null);
        }

        [Test]
        public void ParseFragmentAnnotation_Default_ParsesMultiChargeNeutralLoss()
        {
            var model = new TestableFragmentIntensityModel();
            var result = model.TestParseFragmentAnnotation("y10-H2O+3");

            Assert.That(result.FragmentIdentifier, Is.EqualTo("y10-H2O"));
            Assert.That(result.Charge, Is.EqualTo(3));
            Assert.That(result.NeutralLossFormula, Is.EqualTo("H2O"));
            Assert.That(result.BaseFragmentIdentifier, Is.EqualTo("y10"));
        }

        [Test]
        public void ParseFragmentAnnotation_Altimeter_ExtractsNeutralLossNoCharge()
        {
            var model = new TestableAltimeterModel();
            var result = model.TestParseFragmentAnnotation("IY-NH3");

            Assert.That(result.FragmentIdentifier, Is.EqualTo("IY-NH3"));
            Assert.That(result.Charge, Is.EqualTo(1));
            Assert.That(result.NeutralLossFormula, Is.EqualTo("NH3"));
            Assert.That(result.NeutralLossMass, Is.EqualTo(ChemicalFormula.ParseFormula("NH3").MonoisotopicMass).Within(0.001));
            Assert.That(result.BaseFragmentIdentifier, Is.EqualTo("IY"));
        }

        [Test]
        public void ParseFragmentAnnotation_Altimeter_ExtractsNeutralLossWithCharge()
        {
            var model = new TestableAltimeterModel();
            var result = model.TestParseFragmentAnnotation("IR-NH3^2");

            Assert.That(result.FragmentIdentifier, Is.EqualTo("IR-NH3"));
            Assert.That(result.Charge, Is.EqualTo(2));
            Assert.That(result.NeutralLossFormula, Is.EqualTo("NH3"));
            Assert.That(result.NeutralLossMass, Is.EqualTo(ChemicalFormula.ParseFormula("NH3").MonoisotopicMass).Within(0.001));
            Assert.That(result.BaseFragmentIdentifier, Is.EqualTo("IR"));
        }

        [Test]
        public void ParseFragmentAnnotation_Altimeter_StillThrowsOnNeutralGain()
        {
            var model = new TestableAltimeterModel();
            Assert.Throws<Exception>(() => model.TestParseFragmentAnnotation("IN+CO"));
        }

        [Test]
        public void ParseFragmentAnnotation_UniSpec_ExtractsNeutralLoss()
        {
            var model = new TestableUniSpecModel();
            var result = model.TestParseFragmentAnnotation("y3-H2O+1");

            Assert.That(result.FragmentIdentifier, Is.EqualTo("y3-H2O"));
            Assert.That(result.Charge, Is.EqualTo(1));
            Assert.That(result.NeutralLossFormula, Is.EqualTo("H2O"));
            Assert.That(result.NeutralLossMass, Is.EqualTo(ChemicalFormula.ParseFormula("H2O").MonoisotopicMass).Within(0.001));
            Assert.That(result.BaseFragmentIdentifier, Is.EqualTo("y3"));
        }

        [Test]
        public void ParseFragmentAnnotation_UniSpec_NoChargeNoNeutralLoss()
        {
            var model = new TestableUniSpecModel();
            var result = model.TestParseFragmentAnnotation("a2");

            Assert.That(result.FragmentIdentifier, Is.EqualTo("a2"));
            Assert.That(result.Charge, Is.EqualTo(1));
            Assert.That(result.NeutralLossFormula, Is.Null);
            Assert.That(result.NeutralLossMass, Is.Null);
            Assert.That(result.BaseFragmentIdentifier, Is.EqualTo("a2"));
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // ResponseToPredictions Tests (Intensity Extraction Bug Fix Regression)
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        public void ResponseToPredictions_ExtractsIntensitiesFromCorrectOutput()
        {
            var model = new TestableFragmentIntensityModel(fragmentIonMappingMode: FragmentIonMappingMode.MapToValidatedFullSequence);
            var inputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 35, null, null) { ValidatedFullSequence = "PEPTIDEK" }
            };

            var response = BuildJsonResponse(
                annotations: new[] { "b1+1", "y1+1" },
                mz: new[] { 150.0, 200.0 },
                intensities: new[] { 0.8, 0.6 }
            );

            var predictions = model.TestResponseToPredictions(new[] { response }, inputs);

            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].FragmentIntensities, Is.EquivalentTo(new[] { 0.8, 0.6 }));
            Assert.That(predictions[0].FragmentMZs, Is.EquivalentTo(new[] { 150.0, 200.0 }));
            Assert.That(predictions[0].FragmentAnnotations, Is.EquivalentTo(new[] { "b1+1", "y1+1" }));
        }

        [Test]
        public void ResponseToPredictions_SkipsImpossibleIons()
        {
            var model = new TestableFragmentIntensityModel(fragmentIonMappingMode: FragmentIonMappingMode.MapToValidatedFullSequence);
            var inputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 35, null, null) { ValidatedFullSequence = "PEPTIDEK" }
            };

            var response = BuildJsonResponse(
                annotations: new[] { "b1+1", "y1+1", "a1+1" },
                mz: new[] { 150.0, 200.0, 100.0 },
                intensities: new[] { 0.8, -1.0, 0.6 }
            );

            var predictions = model.TestResponseToPredictions(new[] { response }, inputs);

            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].FragmentIntensities.Count, Is.EqualTo(2));
            Assert.That(predictions[0].FragmentIntensities, Is.EquivalentTo(new[] { 0.8, 0.6 }));
        }

        [Test]
        public void ResponseToPredictions_HandlesMultiplePeptidesInSingleBatch()
        {
            var model = new TestableFragmentIntensityModel(fragmentIonMappingMode: FragmentIonMappingMode.MapToValidatedFullSequence);
            var inputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 35, null, null) { ValidatedFullSequence = "PEPTIDEK" },
                new("ANOTHERSEQ", 2, 35, null, null) { ValidatedFullSequence = "ANOTHERSEQ" }
            };

            var response = BuildJsonResponse(
                annotations: new[] { "b1+1", "y1+1" },
                mz: new[] { 150.0, 200.0 },
                intensities: new[] { 0.8, 0.6 }
            );

            var predictions = model.TestResponseToPredictions(new[] { response }, inputs);

            Assert.That(predictions.Count, Is.EqualTo(2));
            Assert.That(predictions[0].FragmentIntensities, Is.EquivalentTo(new[] { 0.8 }));
            Assert.That(predictions[1].FragmentIntensities, Is.EquivalentTo(new[] { 0.6 }));
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // Helpers
        // ═══════════════════════════════════════════════════════════════════════════

        private static string BuildJsonResponse(string[] annotations, double[] mz, double[] intensities)
        {
            var annotationsJson = string.Join(",", annotations.Select(a => $"\"{a}\""));
            var mzJson = string.Join(",", mz.Select(v => v.ToString(System.Globalization.CultureInfo.InvariantCulture)));
            var intensitiesJson = string.Join(",", intensities.Select(v => v.ToString(System.Globalization.CultureInfo.InvariantCulture)));

            return $"{{\"outputs\":[{{\"name\":\"annotation\",\"datatype\":\"BYTES\",\"shape\":[{annotations.Length}],\"data\":[{annotationsJson}]}},{{\"name\":\"mz\",\"datatype\":\"FP32\",\"shape\":[{mz.Length}],\"data\":[{mzJson}]}},{{\"name\":\"intensities\",\"datatype\":\"FP32\",\"shape\":[{intensities.Length}],\"data\":[{intensitiesJson}]}}]}}";
        }

        /// <summary>
        /// Test harness that exposes protected methods for unit testing.
        /// </summary>
        private sealed class TestableFragmentIntensityModel : FragmentIntensityModel
        {
            private static readonly IReadOnlySet<int> SupportedUnimodIds = new HashSet<int>();
            private static readonly ISequenceConverter Converter = CreateUnimodConverter(
                UnimodSequenceFormatSchema.Instance, SupportedUnimodIds);

            public override string ModelName => "TestModel";
            public override int MaxBatchSize => 1000;
            public override int MaxNumberOfBatchesPerRequest { get; init; } = 1;
            public override int ThrottlingDelayInMilliseconds { get; init; } = 0;
            public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 0;
            public override int MaxPeptideLength => 50;
            public override int MinPeptideLength => 1;
            public override SequenceConversionHandlingMode ModHandlingMode { get; init; }
            public override IncompatibleParameterHandlingMode ParameterHandlingMode { get; init; }
            public override FragmentIonMappingMode FragmentIonMappingMode { get; init; }
            public override HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6 };
            public override HashSet<int> AllowedCollisionEnergies => Enumerable.Range(0, 101).ToHashSet();

            public TestableFragmentIntensityModel(
                FragmentIonMappingMode fragmentIonMappingMode = FragmentIonMappingMode.MapToValidatedFullSequence)
                : base(Converter)
            {
                FragmentIonMappingMode = fragmentIonMappingMode;
            }

            protected override List<Dictionary<string, object>> ToBatchedRequests(List<FragmentIntensityPredictionInput> validInputs)
                => new();

            public (List<object> annotations, List<object> mz, List<object> intensities) TestExtractOutputs(ResponseJSONStruct response)
                => ExtractOutputs(response);

            public ParsedFragmentAnnotation TestParseFragmentAnnotation(string annotation)
                => ParseFragmentAnnotation(annotation);

            public List<PeptideFragmentIntensityPrediction> TestResponseToPredictions(
                IReadOnlyList<string> responses, List<FragmentIntensityPredictionInput> requestInputs)
            {
                ModelInputs = requestInputs;
                return ResponseToPredictions(responses, requestInputs);
            }
        }

        /// <summary>
        /// Test harness for Altimeter2024Intensities that exposes protected methods.
        /// </summary>
        private sealed class TestableAltimeterModel : Altimeter2024Intensities
        {
            public TestableAltimeterModel() : base() { }

            public ParsedFragmentAnnotation TestParseFragmentAnnotation(string annotation)
                => ParseFragmentAnnotation(annotation);
        }

        /// <summary>
        /// Test harness for UniSpec that exposes protected methods.
        /// </summary>
        private sealed class TestableUniSpecModel : UniSpec
        {
            public TestableUniSpecModel() : base() { }

            public ParsedFragmentAnnotation TestParseFragmentAnnotation(string annotation)
                => ParseFragmentAnnotation(annotation);
        }
    }
}
