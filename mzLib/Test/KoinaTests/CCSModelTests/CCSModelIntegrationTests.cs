using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.CCSModels;

namespace Test.KoinaTests.CCSModelTests
{
    [TestFixture]
    [Category("Koina")]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class CCSModelIntegrationTests
    {
        // ═══════════════════════════════════════════════════════════════════════════
        // IM2Deep
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        public static void TestIM2Deep_AcceptsValidPeptides()
        {
            var modelInputs = new List<CCSPredictionInput>
            {
                new("PEPTIDEK", 2),
                new("ACDEFGHIKLMNPQRSTVWY", 3),
                new("LMNPQRSTVWY", 2)
            };

            var model = new IM2Deep();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(3));
            Assert.That(predictions.All(p => p.PredictedCCS.HasValue), Is.True);
        }

        [Test]
        public static void TestIM2Deep_ChargeStateBoundaries()
        {
            var model = new IM2Deep();

            // IM2Deep declares no charge constraint, but positive charges are required.
            for (int charge = 1; charge <= 9; charge++)
            {
                var inputs = new List<CCSPredictionInput> { new("PEPTIDEK", charge) };
                var predictions = model.Predict(inputs);
                Assert.That(predictions[0].PredictedCCS.HasValue, Is.True, $"Charge {charge} should be accepted (no model-level constraint)");
            }

            // Non-physical charges must be rejected client-side, even with no model constraint.
            foreach (int charge in new[] { -1, 0 })
            {
                var inputs = new List<CCSPredictionInput> { new("PEPTIDEK", charge) };
                var predictions = model.Predict(inputs);
                Assert.That(predictions[0].PredictedCCS.HasValue, Is.False, $"Charge {charge} should be rejected as non-physical");
            }
        }

        [Test]
        public static void TestIM2Deep_HandlesSupportedModifications()
        {
            var modelInputs = new List<CCSPredictionInput>
            {
                new("M[Common Variable:Oxidation on M]PEPTIDEK", 2),
                new("C[Common Fixed:Carbamidomethyl on C]PEPTIDEK", 2),
                new("K[Common Variable:Acetyl on K]PEPTIDE", 2)
            };

            var model = new IM2Deep();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(3));
            Assert.That(predictions.All(p => p.PredictedCCS.HasValue), Is.True);
        }

        [Test]
        public static void TestIM2Deep_EmptyInputHandling()
        {
            var model = new IM2Deep();
            var predictions = model.Predict(new List<CCSPredictionInput>());
            Assert.That(predictions.Count, Is.EqualTo(0));
        }

        [Test]
        public static void TestIM2Deep_PredictionQuality()
        {
            var modelInputs = new List<CCSPredictionInput>
            {
                new("PEPTIDEK", 2)
            };

            var model = new IM2Deep();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].PredictedCCS, Is.GreaterThan(0), "CCS value should be positive");
        }

        [Test]
        public static void TestIM2Deep_ModelProperties()
        {
            var model = new IM2Deep();
            Assert.That(model.ModelName, Is.EqualTo("IM2Deep"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(50));
            Assert.That(model.MinPeptideLength, Is.EqualTo(1));
            Assert.That(model.AllowedPrecursorCharges, Is.Empty, "IM2Deep does not constrain precursor charges");
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // AlphaPeptDeepCCSGeneric
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        public static void TestAlphaPeptDeepCCSGeneric_AcceptsValidPeptides()
        {
            var modelInputs = new List<CCSPredictionInput>
            {
                new("PEPTIDEK", 2),
                new("ACDEFGHIKLMNPQRSTVWY", 3),
                new("LMNPQRSTVWY", 2)
            };

            var model = new AlphaPeptDeepCCSGeneric();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(3));
            Assert.That(predictions.All(p => p.PredictedCCS.HasValue), Is.True);
        }

        [Test]
        public static void TestAlphaPeptDeepCCSGeneric_ChargeStateBoundaries()
        {
            var model = new AlphaPeptDeepCCSGeneric();

            for (int charge = 1; charge <= 6; charge++)
            {
                var inputs = new List<CCSPredictionInput> { new("PEPTIDEK", charge) };
                var predictions = model.Predict(inputs);
                Assert.That(predictions[0].PredictedCCS.HasValue, Is.True, $"Charge {charge} should be valid");
            }

            // No charge constraint — charges outside 1-6 should also be accepted
            for (int charge = 7; charge <= 9; charge++)
            {
                var inputs = new List<CCSPredictionInput> { new("PEPTIDEK", charge) };
                var predictions = model.Predict(inputs);
                Assert.That(predictions[0].PredictedCCS.HasValue, Is.True, $"Charge {charge} should be valid (no charge constraint)");
            }

            // Non-physical charges must be rejected client-side, even with no model constraint.
            foreach (int charge in new[] { -1, 0 })
            {
                var inputs = new List<CCSPredictionInput> { new("PEPTIDEK", charge) };
                var predictions = model.Predict(inputs);
                Assert.That(predictions[0].PredictedCCS.HasValue, Is.False, $"Charge {charge} should be rejected as non-physical");
            }
        }

        [Test]
        public static void TestAlphaPeptDeepCCSGeneric_HandlesSupportedModifications()
        {
            var modelInputs = new List<CCSPredictionInput>
            {
                new("M[Common Variable:Oxidation on M]PEPTIDEK", 2),
                new("C[Common Fixed:Carbamidomethyl on C]PEPTIDEK", 2)
            };

            var model = new AlphaPeptDeepCCSGeneric();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(2));
            Assert.That(predictions.All(p => p.PredictedCCS.HasValue), Is.True);
        }

        [Test]
        public static void TestAlphaPeptDeepCCSGeneric_EmptyInputHandling()
        {
            var model = new AlphaPeptDeepCCSGeneric();
            var predictions = model.Predict(new List<CCSPredictionInput>());
            Assert.That(predictions.Count, Is.EqualTo(0));
        }

        [Test]
        public static void TestAlphaPeptDeepCCSGeneric_PredictionQuality()
        {
            var modelInputs = new List<CCSPredictionInput>
            {
                new("PEPTIDEK", 2)
            };

            var model = new AlphaPeptDeepCCSGeneric();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].PredictedCCS, Is.GreaterThan(0), "CCS value should be positive");
        }

        [Test]
        public static void TestAlphaPeptDeepCCSGeneric_ModelProperties()
        {
            var model = new AlphaPeptDeepCCSGeneric();
            Assert.That(model.ModelName, Is.EqualTo("AlphaPeptDeep_ccs_generic"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(50));
            Assert.That(model.MinPeptideLength, Is.EqualTo(1));
            Assert.That(model.AllowedPrecursorCharges, Is.Empty, "AlphaPeptDeepCCSGeneric does not constrain precursor charges");
        }
    }
}
