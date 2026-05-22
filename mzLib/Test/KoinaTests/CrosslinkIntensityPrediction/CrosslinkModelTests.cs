using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.CrosslinkIntensityModels;
using PredictionClients.Koina.Util;

namespace Test.KoinaTests.CrosslinkIntensityPrediction
{
    [TestFixture]
    [Category("Koina")]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class CrosslinkModelTests
    {
        // ═══════════════════════════════════════════════════════════════════════════
        // Prosit2023IntensityXLCMS2
        // ═══════════════════════════════════════════════════════════════════════════

        private const string Cms2AlphaSeq = "PEPTIDEK[UNIMOD:1896]";
        private const string Cms2BetaSeq = "ACDEFGK[UNIMOD:1896]HIK";
        private const string Cms2AlphaSeq2 = "LMNPEK[UNIMOD:1896]QRSTVWY";
        private const string Cms2BetaSeq2 = "PEPK[UNIMOD:1896]TIDEK";
        private const string Cms3AlphaSeq = "PEPTIDEK[UNIMOD:1881]";
        private const string Cms3BetaSeq = "ACDEFGK[UNIMOD:1881]HIK";
        private const string Cms3AlphaSeq2 = "LMNPEK[UNIMOD:1881]QRSTVWY";
        private const string Cms3BetaSeq2 = "PEPK[UNIMOD:1881]TIDEK";
        private const string Nms2AlphaSeq = "PEPTIDEK[UNIMOD:1898]";
        private const string Nms2BetaSeq = "ACDEFGK[UNIMOD:1898]HIK";
        private const string Nms2AlphaSeq2 = "LMNPEK[UNIMOD:1898]QRSTVWY";
        private const string Nms2BetaSeq2 = "PEPK[UNIMOD:1898]TIDEK";

        [Test]
        public static void TestProsit2023XLCMS2_AcceptsValidPeptides()
        {
            var modelInputs = new List<CrosslinkIntensityPredictionInput>
            {
                new(Cms2AlphaSeq, Cms2BetaSeq, 2, 35),
                new(Cms2AlphaSeq2, Cms2BetaSeq2, 3, 35)
            };

            var model = new Prosit2023IntensityXLCMS2();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(2));
            Assert.That(predictions.All(p => p.FragmentAnnotations != null), Is.True);
        }

        [Test]
        public static void TestProsit2023XLCMS2_ChargeStateBoundaries()
        {
            var model = new Prosit2023IntensityXLCMS2();

            for (int charge = 1; charge <= 6; charge++)
            {
                var inputs = new List<CrosslinkIntensityPredictionInput>
                {
                    new(Cms2AlphaSeq, Cms2BetaSeq, charge, 35)
                };
                var predictions = model.Predict(inputs);
                Assert.That(predictions[0].FragmentAnnotations, Is.Not.Null, $"Charge {charge} should be valid");
            }

            var invalidInputs = new List<CrosslinkIntensityPredictionInput>
            {
                new(Cms2AlphaSeq, Cms2BetaSeq, 0, 35)
            };
            var invalidPredictions = model.Predict(invalidInputs);
            Assert.That(invalidPredictions[0].FragmentAnnotations, Is.Null, "Charge 0 should be invalid");
        }

        [Test]
        public static void TestProsit2023XLCMS2_EmptyInputHandling()
        {
            var model = new Prosit2023IntensityXLCMS2();
            var predictions = model.Predict(new List<CrosslinkIntensityPredictionInput>());
            Assert.That(predictions.Count, Is.EqualTo(0));
        }

        [Test]
        public static void TestProsit2023XLCMS2_PredictionQuality()
        {
            var modelInputs = new List<CrosslinkIntensityPredictionInput>
            {
                new(Cms2AlphaSeq, Cms2BetaSeq, 2, 35)
            };

            var model = new Prosit2023IntensityXLCMS2();
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
        public static void TestProsit2023XLCMS2_ModelProperties()
        {
            var model = new Prosit2023IntensityXLCMS2();
            Assert.That(model.ModelName, Is.EqualTo("Prosit_2023_intensity_XL_CMS2"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(30));
            Assert.That(model.NumberOfPredictedFragmentIons, Is.EqualTo(348));
            Assert.That(model.AllowedPrecursorCharges, Is.EquivalentTo(new[] { 1, 2, 3, 4, 5, 6 }));
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // Prosit2023IntensityXLCMS3
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        public static void TestProsit2023XLCMS3_AcceptsValidPeptides()
        {
            var modelInputs = new List<CrosslinkIntensityPredictionInput>
            {
                new(Cms3AlphaSeq, Cms3BetaSeq, 2, 35),
                new(Cms3AlphaSeq2, Cms3BetaSeq2, 3, 35)
            };

            var model = new Prosit2023IntensityXLCMS3();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(2));
            Assert.That(predictions.All(p => p.FragmentAnnotations != null), Is.True);
        }

        [Test]
        public static void TestProsit2023XLCMS3_ChargeStateBoundaries()
        {
            var model = new Prosit2023IntensityXLCMS3();

            for (int charge = 1; charge <= 6; charge++)
            {
                var inputs = new List<CrosslinkIntensityPredictionInput>
                {
                    new(Cms3AlphaSeq, Cms3BetaSeq, charge, 35)
                };
                var predictions = model.Predict(inputs);
                Assert.That(predictions[0].FragmentAnnotations, Is.Not.Null, $"Charge {charge} should be valid");
            }
        }

        [Test]
        public static void TestProsit2023XLCMS3_EmptyInputHandling()
        {
            var model = new Prosit2023IntensityXLCMS3();
            var predictions = model.Predict(new List<CrosslinkIntensityPredictionInput>());
            Assert.That(predictions.Count, Is.EqualTo(0));
        }

        [Test]
        public static void TestProsit2023XLCMS3_PredictionQuality()
        {
            var modelInputs = new List<CrosslinkIntensityPredictionInput>
            {
                new(Cms3AlphaSeq, Cms3BetaSeq, 2, 35)
            };

            var model = new Prosit2023IntensityXLCMS3();
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
        public static void TestProsit2023XLCMS3_ModelProperties()
        {
            var model = new Prosit2023IntensityXLCMS3();
            Assert.That(model.ModelName, Is.EqualTo("Prosit_2023_intensity_XL_CMS3"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(30));
            Assert.That(model.NumberOfPredictedFragmentIons, Is.EqualTo(174));
            Assert.That(model.AllowedCollisionEnergies.Count, Is.EqualTo(0), "CMS3 uses fixed NCE");
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // Prosit2024IntensityXLNMS2
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        public static void TestProsit2024XLNMS2_AcceptsValidPeptides()
        {
            var modelInputs = new List<CrosslinkIntensityPredictionInput>
            {
                new(Nms2AlphaSeq, Nms2BetaSeq, 2, 35),
                new(Nms2AlphaSeq2, Nms2BetaSeq2, 3, 35)
            };

            var model = new Prosit2024IntensityXLNMS2();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(2));
            Assert.That(predictions.All(p => p.FragmentAnnotations != null), Is.True);
        }

        [Test]
        public static void TestProsit2024XLNMS2_ChargeStateBoundaries()
        {
            var model = new Prosit2024IntensityXLNMS2();

            for (int charge = 1; charge <= 6; charge++)
            {
                var inputs = new List<CrosslinkIntensityPredictionInput>
                {
                    new(Nms2AlphaSeq, Nms2BetaSeq, charge, 35)
                };
                var predictions = model.Predict(inputs);
                Assert.That(predictions[0].FragmentAnnotations, Is.Not.Null, $"Charge {charge} should be valid");
            }
        }

        [Test]
        public static void TestProsit2024XLNMS2_EmptyInputHandling()
        {
            var model = new Prosit2024IntensityXLNMS2();
            var predictions = model.Predict(new List<CrosslinkIntensityPredictionInput>());
            Assert.That(predictions.Count, Is.EqualTo(0));
        }

        [Test]
        public static void TestProsit2024XLNMS2_PredictionQuality()
        {
            var modelInputs = new List<CrosslinkIntensityPredictionInput>
            {
                new(Nms2AlphaSeq, Nms2BetaSeq, 2, 35)
            };

            var model = new Prosit2024IntensityXLNMS2();
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
        public static void TestProsit2024XLNMS2_AllowedCollisionEnergies()
        {
            var model = new Prosit2024IntensityXLNMS2();

            // Koina CMS2 accepts any CE, so empty set means no restriction
            Assert.That(model.AllowedCollisionEnergies, Is.Empty);
        }

        [Test]
        public static void TestProsit2024XLNMS2_ModelProperties()
        {
            var model = new Prosit2024IntensityXLNMS2();
            Assert.That(model.ModelName, Is.EqualTo("Prosit_2024_intensity_XL_NMS2"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(30));
            Assert.That(model.NumberOfPredictedFragmentIons, Is.EqualTo(174));
            Assert.That(model.AllowedPrecursorCharges, Is.EquivalentTo(new[] { 1, 2, 3, 4, 5, 6 }));
            Assert.That(model.AllowedCollisionEnergies.Count, Is.EqualTo(0), "NMS2 accepts any FP32 collision energy value");
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // Modification & sequence validation (no network: all inputs rejected client-side)
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        public static void TestCrosslinkModel_RejectsDisallowedUnimodId()
        {
            // UNIMOD:21 (phospho) is not in any crosslink model's allowed set.
            var model = new Prosit2023IntensityXLCMS2();
            var inputs = new List<CrosslinkIntensityPredictionInput>
            {
                new("PEPS[UNIMOD:21]IDEK[UNIMOD:1896]", "ACDEFGK[UNIMOD:1896]HIK", 2, 35)
            };
            var predictions = model.Predict(inputs);
            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].FragmentAnnotations, Is.Null);
            Assert.That(predictions[0].Warning, Is.Not.Null);
            Assert.That(predictions[0].Warning!.Message, Does.Contain("UNIMOD:21").Or.Contain("21"));
        }

        [Test]
        public static void TestCrosslinkModel_RejectsNonUnimodNotation()
        {
            // The Koina helper expects K[UNIMOD:1896], not mzLib's "Common Variable:..." notation.
            var model = new Prosit2023IntensityXLCMS2();
            var inputs = new List<CrosslinkIntensityPredictionInput>
            {
                new("K[Common Variable:Acetyl on K]PEPTIDEK[UNIMOD:1896]", "ACDEFGK[UNIMOD:1896]HIK", 2, 35)
            };
            var predictions = model.Predict(inputs);
            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].FragmentAnnotations, Is.Null);
            Assert.That(predictions[0].Warning, Is.Not.Null);
        }

        [Test]
        public static void TestCrosslinkModel_ThrowOnDisallowedUnimodId()
        {
            var model = new Prosit2023IntensityXLCMS2(
                modHandlingMode: SequenceConversionHandlingMode.ThrowException);
            var inputs = new List<CrosslinkIntensityPredictionInput>
            {
                new("PEPS[UNIMOD:21]IDEK[UNIMOD:1896]", "ACDEFGK[UNIMOD:1896]HIK", 2, 35)
            };
            Assert.Throws<ArgumentException>(() => model.Predict(inputs));
        }

        [Test]
        public static void TestPairedCrosslinkModel_RejectsNullBeta()
        {
            // CMS2 requires both alpha and beta.
            var model = new Prosit2023IntensityXLCMS2();
            var inputs = new List<CrosslinkIntensityPredictionInput>
            {
                new(Cms2AlphaSeq, null, 2, 35)
            };
            var predictions = model.Predict(inputs);
            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].FragmentAnnotations, Is.Null);
            Assert.That(predictions[0].Warning, Is.Not.Null);
            Assert.That(predictions[0].Warning!.Message, Does.Contain("Beta sequence is required").IgnoreCase);
        }

        [Test]
        public static void TestSingleSequenceCrosslinkModel_AcceptsNullBeta()
        {
            // CMS3 takes a single combined sequence; beta is unused and may be null.
            // Validate this without hitting Koina by inspecting the model property.
            var model = new Prosit2023IntensityXLCMS3();
            Assert.That(model.RequiresBetaSequence, Is.False);

            // Verify other crosslink models still require beta.
            Assert.That(new Prosit2023IntensityXLCMS2().RequiresBetaSequence, Is.True);
            Assert.That(new Prosit2024IntensityXLNMS2().RequiresBetaSequence, Is.True);
        }
    }
}
