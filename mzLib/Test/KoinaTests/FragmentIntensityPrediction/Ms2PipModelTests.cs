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
    public class Ms2PipModelTests
    {
        // ═══════════════════════════════════════════════════════════════════════════
        // Shared Tests for All MS2PIP Models
        // ═══════════════════════════════════════════════════════════════════════════

        private static void TestMs2PipModelAcceptsValidPeptides(FragmentIntensityModel model)
        {
            var modelInputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 35, null, null),
                new("ACDEFGHIKLMNPQRSTVWY", 3, 35, null, null),
                new("LMNPQRSTVWY", 2, 35, null, null)
            };

            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(3), $"Model {model.ModelName} should return predictions for all inputs");
            Assert.That(predictions.All(p => p.FragmentAnnotations != null), Is.True, $"Model {model.ModelName} should have predictions for all valid peptides");
        }

        private static void TestMs2PipModelChargeStates(FragmentIntensityModel model)
        {
            var throwModel = CreateThrowModel(model);

            for (int charge = 1; charge <= 6; charge++)
            {
                var inputs = new List<FragmentIntensityPredictionInput>
                {
                    new("PEPTIDEK", charge, 35, null, null)
                };
                Assert.DoesNotThrow(() => throwModel.Predict(inputs), $"Model {model.ModelName}: Charge {charge} should be valid");
            }

            var invalidInputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 0, 35, null, null)
            };
            Assert.Throws<ArgumentException>(() => throwModel.Predict(invalidInputs), $"Model {model.ModelName}: Charge 0 should throw");
        }

        private static void TestMs2PipModelPredictionQuality(FragmentIntensityModel model)
        {
            var modelInputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 35, null, null)
            };

            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(1), $"Model {model.ModelName} should return prediction");
            var prediction = predictions[0];

            Assert.That(prediction.FragmentAnnotations!.Count, Is.GreaterThan(0), $"Model {model.ModelName} should have fragment ions");
            for (int i = 0; i < prediction.FragmentAnnotations.Count; i++)
            {
                Assert.That(prediction.FragmentMZs![i], Is.GreaterThan(0), $"Model {model.ModelName}: m/z should be positive");
                Assert.That(prediction.FragmentIntensities![i], Is.GreaterThanOrEqualTo(-1), $"Model {model.ModelName}: intensity should be >= -1");
            }
        }

        private static void TestMs2PipModelEmptyInput(FragmentIntensityModel model)
        {
            var predictions = model.Predict(new List<FragmentIntensityPredictionInput>());
            Assert.That(predictions.Count, Is.EqualTo(0), $"Model {model.ModelName}: empty input should return empty predictions");
        }

        private static FragmentIntensityModel CreateThrowModel(FragmentIntensityModel model)
        {
            return model switch
            {
                Ms2PipHCD2021 m => new Ms2PipHCD2021(parameterHandlingMode: IncompatibleParameterHandlingMode.ThrowException),
                Ms2PipCIDTMT m => new Ms2PipCIDTMT(parameterHandlingMode: IncompatibleParameterHandlingMode.ThrowException),
                Ms2PipImmunoHCD m => new Ms2PipImmunoHCD(parameterHandlingMode: IncompatibleParameterHandlingMode.ThrowException),
                Ms2PipTTOF5600 m => new Ms2PipTTOF5600(parameterHandlingMode: IncompatibleParameterHandlingMode.ThrowException),
                Ms2PipITRAQPhospho m => new Ms2PipITRAQPhospho(parameterHandlingMode: IncompatibleParameterHandlingMode.ThrowException),
                Ms2PipTimsTOF2023 m => new Ms2PipTimsTOF2023(parameterHandlingMode: IncompatibleParameterHandlingMode.ThrowException),
                Ms2PipTimsTOF2024 m => new Ms2PipTimsTOF2024(parameterHandlingMode: IncompatibleParameterHandlingMode.ThrowException),
                _ => throw new ArgumentException($"Unknown model type: {model.GetType().Name}")
            };
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // ms2pip_HCD2021
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        public static void TestMs2PipHCD2021_AcceptsValidPeptides()
        {
            TestMs2PipModelAcceptsValidPeptides(new Ms2PipHCD2021());
        }

        [Test]
        public static void TestMs2PipHCD2021_ChargeStates()
        {
            TestMs2PipModelChargeStates(new Ms2PipHCD2021());
        }

        [Test]
        public static void TestMs2PipHCD2021_PredictionQuality()
        {
            TestMs2PipModelPredictionQuality(new Ms2PipHCD2021());
        }

        [Test]
        public static void TestMs2PipHCD2021_EmptyInput()
        {
            TestMs2PipModelEmptyInput(new Ms2PipHCD2021());
        }

        [Test]
        public static void TestMs2PipHCD2021_ModelProperties()
        {
            var model = new Ms2PipHCD2021();
            Assert.That(model.ModelName, Is.EqualTo("ms2pip_HCD2021"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(30));
            Assert.That(model.NumberOfPredictedFragmentIons, Is.EqualTo(58));
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // ms2pip_CID_TMT
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        public static void TestMs2PipCIDTMT_AcceptsValidPeptides()
        {
            var modelInputs = new List<FragmentIntensityPredictionInput>
            {
                new("[Common Fixed:TMT6plex on N-terminus]PEPTIDEK", 2, 35, null, null),
                new("[Common Fixed:TMT6plex on N-terminus]ACDEFGHIKLMNPQRSTVWY", 3, 35, null, null)
            };

            var model = new Ms2PipCIDTMT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(2));
            Assert.That(predictions.All(p => p.FragmentAnnotations != null), Is.True);
        }

        [Test]
        public static void TestMs2PipCIDTMT_ChargeStates()
        {
            TestMs2PipModelChargeStates(new Ms2PipCIDTMT());
        }

        [Test]
        public static void TestMs2PipCIDTMT_PredictionQuality()
        {
            var modelInputs = new List<FragmentIntensityPredictionInput>
            {
                new("[Common Fixed:TMT6plex on N-terminus]PEPTIDEK", 2, 35, null, null)
            };

            var model = new Ms2PipCIDTMT();
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

        [Test]
        public static void TestMs2PipCIDTMT_ModelProperties()
        {
            var model = new Ms2PipCIDTMT();
            Assert.That(model.ModelName, Is.EqualTo("ms2pip_CID_TMT"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(30));
            Assert.That(model.NumberOfPredictedFragmentIons, Is.EqualTo(58));
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // ms2pip_Immuno_HCD
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        public static void TestMs2PipImmunoHCD_AcceptsValidPeptides()
        {
            var modelInputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 35, null, null),
                new("ACDEFGHIK", 3, 35, null, null)
            };

            var model = new Ms2PipImmunoHCD();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(2));
            Assert.That(predictions.All(p => p.FragmentAnnotations != null), Is.True);
        }

        [Test]
        public static void TestMs2PipImmunoHCD_PeptideLengthBoundary()
        {
            var model = new Ms2PipImmunoHCD();
            Assert.That(model.MaxPeptideLength, Is.EqualTo(30));

            // Well under the boundary (13 AA).
            var validInput = new List<FragmentIntensityPredictionInput>
            {
                new("ACDEFGHIKLMNP", 2, 35, null, null)
            };
            var predictions = model.Predict(validInput);
            Assert.That(predictions.Count, Is.EqualTo(1));

            // Exactly at the boundary (30 AA) — must be accepted.
            var atBoundary = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKL"; // 30 chars
            Assert.That(atBoundary.Length, Is.EqualTo(30));
            var boundaryInput = new List<FragmentIntensityPredictionInput>
            {
                new(atBoundary, 2, 35, null, null)
            };
            predictions = model.Predict(boundaryInput);
            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].FragmentAnnotations, Is.Not.Null, "Peptide at the 30 AA boundary should be accepted");

            // One over the boundary (31 AA) — must be filtered.
            var overBoundary = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLM"; // 31 chars
            Assert.That(overBoundary.Length, Is.EqualTo(31));
            var invalidInput = new List<FragmentIntensityPredictionInput>
            {
                new(overBoundary, 2, 35, null, null)
            };
            predictions = model.Predict(invalidInput);
            Assert.That(predictions[0].FragmentAnnotations, Is.Null, "Peptide > 30 AA should be filtered");
        }

        [Test]
        public static void TestMs2PipImmunoHCD_ModelProperties()
        {
            var model = new Ms2PipImmunoHCD();
            Assert.That(model.ModelName, Is.EqualTo("ms2pip_Immuno_HCD"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(30));
            Assert.That(model.NumberOfPredictedFragmentIons, Is.EqualTo(58));
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // ms2pip_TTOF5600
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        public static void TestMs2PipTTOF5600_AcceptsValidPeptides()
        {
            TestMs2PipModelAcceptsValidPeptides(new Ms2PipTTOF5600());
        }

        [Test]
        public static void TestMs2PipTTOF5600_ChargeStates()
        {
            TestMs2PipModelChargeStates(new Ms2PipTTOF5600());
        }

        [Test]
        public static void TestMs2PipTTOF5600_PredictionQuality()
        {
            TestMs2PipModelPredictionQuality(new Ms2PipTTOF5600());
        }

        [Test]
        public static void TestMs2PipTTOF5600_ModelProperties()
        {
            var model = new Ms2PipTTOF5600();
            Assert.That(model.ModelName, Is.EqualTo("ms2pip_TTOF5600"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(30));
            Assert.That(model.NumberOfPredictedFragmentIons, Is.EqualTo(58));
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // ms2pip_iTRAQphospho
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        public static void TestMs2PipITRAQPhospho_AcceptsValidPeptides()
        {
            TestMs2PipModelAcceptsValidPeptides(new Ms2PipITRAQPhospho());
        }

        [Test]
        public static void TestMs2PipITRAQPhospho_ChargeStates()
        {
            TestMs2PipModelChargeStates(new Ms2PipITRAQPhospho());
        }

        [Test]
        public static void TestMs2PipITRAQPhospho_PredictionQuality()
        {
            TestMs2PipModelPredictionQuality(new Ms2PipITRAQPhospho());
        }

        [Test]
        public static void TestMs2PipITRAQPhospho_ModelProperties()
        {
            var model = new Ms2PipITRAQPhospho();
            Assert.That(model.ModelName, Is.EqualTo("ms2pip_iTRAQphospho"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(30));
            Assert.That(model.NumberOfPredictedFragmentIons, Is.EqualTo(58));
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // ms2pip_timsTOF2023
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        public static void TestMs2PipTimsTOF2023_AcceptsValidPeptides()
        {
            TestMs2PipModelAcceptsValidPeptides(new Ms2PipTimsTOF2023());
        }

        [Test]
        public static void TestMs2PipTimsTOF2023_ChargeStates()
        {
            TestMs2PipModelChargeStates(new Ms2PipTimsTOF2023());
        }

        [Test]
        public static void TestMs2PipTimsTOF2023_PredictionQuality()
        {
            TestMs2PipModelPredictionQuality(new Ms2PipTimsTOF2023());
        }

        [Test]
        public static void TestMs2PipTimsTOF2023_ModelProperties()
        {
            var model = new Ms2PipTimsTOF2023();
            Assert.That(model.ModelName, Is.EqualTo("ms2pip_timsTOF2023"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(30));
            Assert.That(model.NumberOfPredictedFragmentIons, Is.EqualTo(58));
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // ms2pip_timsTOF2024
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        public static void TestMs2PipTimsTOF2024_AcceptsValidPeptides()
        {
            TestMs2PipModelAcceptsValidPeptides(new Ms2PipTimsTOF2024());
        }

        [Test]
        public static void TestMs2PipTimsTOF2024_ChargeStates()
        {
            TestMs2PipModelChargeStates(new Ms2PipTimsTOF2024());
        }

        [Test]
        public static void TestMs2PipTimsTOF2024_PredictionQuality()
        {
            TestMs2PipModelPredictionQuality(new Ms2PipTimsTOF2024());
        }

        [Test]
        public static void TestMs2PipTimsTOF2024_ModelProperties()
        {
            var model = new Ms2PipTimsTOF2024();
            Assert.That(model.ModelName, Is.EqualTo("ms2pip_timsTOF2024"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(30));
            Assert.That(model.NumberOfPredictedFragmentIons, Is.EqualTo(58));
        }
    }
}
