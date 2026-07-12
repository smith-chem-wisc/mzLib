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
    public class Altimeter2024IntensitiesTests
    {
        /// <summary>
        /// Tests that the model correctly processes valid peptide sequences.
        /// </summary>
        [Test]
        public static void TestAltimeterAcceptsValidPeptides()
        {
            var modelInputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 35, null, null),
                new("C[Common Fixed:Carbamidomethyl on C]DEFGHIKLMNPQRSTVWY", 3, 35, null, null),
                new("LMNPQRSTVWY", 2, 35, null, null)
            };

            var model = new Altimeter2024Intensities();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(3));
            Assert.That(predictions.All(p => p.FragmentAnnotations != null), Is.True);
        }

        /// <summary>
        /// Tests that the model correctly handles precursor charge constraints (1-7).
        /// </summary>
        [Test]
        public static void TestAltimeterChargeStateBoundaries()
        {
            var throwModel = new Altimeter2024Intensities(parameterHandlingMode: IncompatibleParameterHandlingMode.ThrowException);

            for (int charge = 1; charge <= 7; charge++)
            {
                var inputs = new List<FragmentIntensityPredictionInput>
                {
                    new("PEPTIDEK", charge, 35, null, null)
                };
                Assert.DoesNotThrow(() => throwModel.Predict(inputs), $"Charge {charge} should be valid");
            }

            var invalidInputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 0, 35, null, null),
                new("PEPTIDEK", 8, 35, null, null)
            };
            Assert.Throws<ArgumentException>(() => throwModel.Predict(invalidInputs));
        }

        /// <summary>
        /// Tests that the model handles modified peptides with supported Unimod IDs (35, 4).
        /// </summary>
        [Test]
        public static void TestAltimeterHandlesSupportedModifications()
        {
            var modelInputs = new List<FragmentIntensityPredictionInput>
            {
                new("M[Common Variable:Oxidation on M]PEPTIDEK", 2, 35, null, null),
                new("C[Common Fixed:Carbamidomethyl on C]PEPTIDEK", 2, 35, null, null)
            };

            var model = new Altimeter2024Intensities();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(2));
            Assert.That(predictions.All(p => p.FragmentAnnotations != null), Is.True);
        }

        /// <summary>
        /// Tests that unsupported modifications are handled according to ModHandlingMode.
        /// </summary>
        [Test]
        public static void TestAltimeterUnsupportedModifications()
        {
            var modelInputs = new List<FragmentIntensityPredictionInput>
            {
                new("S[Common Variable:Phospho on S]PEPTIDEK", 2, 35, null, null)
            };

            var model = new Altimeter2024Intensities(modHandlingMode: SequenceConversionHandlingMode.ReturnNull);
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].FragmentAnnotations, Is.Null, "Unsupported mod (Phospho) should result in null prediction");
        }

        /// <summary>
        /// Tests empty input handling.
        /// </summary>
        [Test]
        public static void TestAltimeterEmptyInputHandling()
        {
            var model = new Altimeter2024Intensities();
            var predictions = model.Predict(new List<FragmentIntensityPredictionInput>());

            Assert.That(predictions.Count, Is.EqualTo(0));
        }

        /// <summary>
        /// Tests that predictions contain valid m/z and intensity values.
        /// </summary>
        [Test]
        public static void TestAltimeterPredictionQuality()
        {
            var modelInputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 35, null, null)
            };

            var model = new Altimeter2024Intensities();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(1));
            var prediction = predictions[0];

            Assert.That(prediction.FragmentAnnotations!.Count, Is.GreaterThan(0));
            for (int i = 0; i < prediction.FragmentAnnotations.Count; i++)
            {
                Assert.That(prediction.FragmentIntensities![i], Is.GreaterThanOrEqualTo(-1));
                if (prediction.FragmentIntensities[i] >= 0)
                {
                    Assert.That(prediction.FragmentMZs![i], Is.GreaterThan(0));
                }
            }
        }

        /// <summary>
        /// Tests model property constants.
        /// </summary>
        [Test]
        public static void TestAltimeterModelProperties()
        {
            var model = new Altimeter2024Intensities();

            Assert.That(model.ModelName, Is.EqualTo("Altimeter_2024_intensities"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(40)); // Koina docs: 6-40
            Assert.That(model.MinPeptideLength, Is.EqualTo(6));
            Assert.That(model.NumberOfPredictedFragmentIons, Is.EqualTo(380));
            Assert.That(model.AllowedPrecursorCharges, Is.EquivalentTo(new[] { 1, 2, 3, 4, 5, 6, 7 }));
        }
    }
}
