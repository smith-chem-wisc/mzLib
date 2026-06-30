using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.Util;

namespace Test.KoinaTests.FragmentIntensityPrediction
{
    [TestFixture]
    [Category("Koina")]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class UniSpecTests
    {
        /// <summary>
        /// Tests that the model correctly processes valid peptide sequences with required instrument type.
        /// </summary>
        [Test]
        public static void TestUniSpecAcceptsValidPeptides()
        {
            var modelInputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 35, "QE", null),
                new("ACDEFGHIKLMNPQRSTVWY", 3, 35, "LUMOS", null),
                new("LMNPQRSTVWY", 2, 35, "NONE", null)
            };

            var model = new UniSpec();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(3));
            Assert.That(predictions.All(p => p.FragmentAnnotations != null), Is.True);
        }

        /// <summary>
        /// Tests that instrument type is required and validated.
        /// </summary>
        [Test]
        public static void TestUniSpecRequiresInstrumentType()
        {
            var modelInputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 35, null, null)
            };

            var model = new UniSpec(parameterHandlingMode: IncompatibleParameterHandlingMode.ReturnNull);
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].FragmentAnnotations, Is.Null, "Missing instrument type should result in null prediction");
        }

        /// <summary>
        /// Tests that unsupported instrument types are rejected.
        /// </summary>
        [Test]
        public static void TestUniSpecRejectsInvalidInstrumentType()
        {
            var modelInputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 35, "INVALID", null)
            };

            var model = new UniSpec(parameterHandlingMode: IncompatibleParameterHandlingMode.ThrowException);
            Assert.Throws<ArgumentException>(() => model.Predict(modelInputs));
        }

        /// <summary>
        /// Tests that the model handles precursor charge constraints with instrument-dependent rules.
        /// QE/QEHFX/ELITE: 2-5, VELOS: 2-4, NONE/LUMOS: any (1-8).
        /// </summary>
        [Test]
        public static void TestUniSpecChargeStateBoundaries()
        {
            var throwModel = new UniSpec(parameterHandlingMode: IncompatibleParameterHandlingMode.ThrowException);

            // NONE allows any charge 1-8
            for (int charge = 1; charge <= 8; charge++)
            {
                var inputs = new List<FragmentIntensityPredictionInput>
                {
                    new("PEPTIDEK", charge, 35, "NONE", null)
                };
                Assert.DoesNotThrow(() => throwModel.Predict(inputs), $"Charge {charge} with NONE should be valid");
            }

            // QE only allows 2-5
            var validQeInputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 35, "QE", null),
                new("PEPTIDEK", 5, 35, "QE", null)
            };
            Assert.DoesNotThrow(() => throwModel.Predict(validQeInputs));

            var invalidQeInputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 1, 35, "QE", null),
                new("PEPTIDEK", 6, 35, "QE", null)
            };
            Assert.Throws<ArgumentException>(() => throwModel.Predict(invalidQeInputs));

            // VELOS only allows 2-4
            var invalidVelosInputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 1, 35, "VELOS", null),
                new("PEPTIDEK", 5, 35, "VELOS", null)
            };
            Assert.Throws<ArgumentException>(() => throwModel.Predict(invalidVelosInputs));
        }

        /// <summary>
        /// Tests that the model handles modified peptides with supported Unimod IDs (35, 4).
        /// </summary>
        [Test]
        public static void TestUniSpecHandlesSupportedModifications()
        {
            var modelInputs = new List<FragmentIntensityPredictionInput>
            {
                new("M[Common Variable:Oxidation on M]PEPTIDEK", 2, 35, "QE", null),
                new("C[Common Fixed:Carbamidomethyl on C]PEPTIDEK", 2, 35, "QE", null)
            };

            var model = new UniSpec();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(2));
            Assert.That(predictions.All(p => p.FragmentAnnotations != null), Is.True);
        }

        /// <summary>
        /// Tests empty input handling.
        /// </summary>
        [Test]
        public static void TestUniSpecEmptyInputHandling()
        {
            var model = new UniSpec();
            var predictions = model.Predict(new List<FragmentIntensityPredictionInput>());

            Assert.That(predictions.Count, Is.EqualTo(0));
        }

        /// <summary>
        /// Tests that predictions contain valid m/z and intensity values.
        /// </summary>
        [Test]
        public static void TestUniSpecPredictionQuality()
        {
            var modelInputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 35, "QE", null)
            };

            var model = new UniSpec();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(1));
            var prediction = predictions[0];

            Assert.That(prediction.FragmentAnnotations!.Count, Is.GreaterThan(0));
            for (int i = 0; i < prediction.FragmentAnnotations.Count; i++)
            {
                Assert.That(prediction.FragmentMZs![i], Is.GreaterThan(0));
                Assert.That(prediction.FragmentIntensities![i], Is.GreaterThanOrEqualTo(-1));
            }
        }

        /// <summary>
        /// Tests model property constants.
        /// </summary>
        [Test]
        public static void TestUniSpecModelProperties()
        {
            var model = new UniSpec();

            Assert.That(model.ModelName, Is.EqualTo("UniSpec"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(40));
            Assert.That(model.MinPeptideLength, Is.EqualTo(1));
            Assert.That(model.AllowedPrecursorCharges, Is.EquivalentTo(new[] { 1, 2, 3, 4, 5, 6, 7, 8 }));
            Assert.That(model.AllowedInstrumentTypes, Is.EquivalentTo(new[] { "QE", "QEHFX", "LUMOS", "ELITE", "VELOS", "NONE" }));
        }
    }
}
