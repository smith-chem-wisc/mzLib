using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using Chemistry;
using NUnit.Framework;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.Util;

namespace Test.KoinaTests.FragmentIntensityPrediction
{
    /// <summary>
    /// No-network coverage for FragmentIntensityModel paths that surround the HTTP call:
    /// the all-invalid Predict realignment, MapToInputFullSequence response mapping, and
    /// spectral-library generation from seeded predictions.
    /// </summary>
    [TestFixture]
    public class FragmentIntensityModelCoverageTests
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

        // ── Predict realignment without HTTP (every input fails validation) ─────────

        [Test]
        public void Predict_AllInputsInvalid_ReturnsPlaceholdersWithWarnings()
        {
            var model = new CoverageModel();
            var inputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEP*TIDE", 2, 30, null, null),   // invalid character
                new("PEPTIDEK", 99, 30, null, null)   // charge out of allowed range
            };

            var predictions = model.Predict(inputs);

            Assert.That(predictions.Count, Is.EqualTo(2));
            Assert.That(predictions.All(p => p.FragmentAnnotations is null), Is.True);
            Assert.That(predictions.All(p => p.Warning is not null), Is.True);
            Assert.That(model.ValidInputsMask, Is.All.False);
            // The invalid input's original sequence/charge are still echoed back.
            Assert.That(predictions[1].PrecursorCharge, Is.EqualTo(99));
        }

        [Test]
        public void Predict_EmptyInput_ReturnsEmpty()
        {
            var model = new CoverageModel();
            Assert.That(model.Predict(new List<FragmentIntensityPredictionInput>()), Is.Empty);
        }

        // ── ResponseToPredictions with MapToInputFullSequence ───────────────────────

        [Test]
        public void ResponseToPredictions_MapToInputFullSequence_RecomputesMzFromInputSequence()
        {
            var model = new CoverageModel(FragmentIonMappingMode.MapToInputFullSequence);
            var inputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 30, null, null) { ValidatedFullSequence = "PEPTIDEK" }
            };
            var response = BuildJsonResponse(
                annotations: new[] { "b2+1", "y3+1" },
                mz: new[] { 0.0, 0.0 },           // ignored: recomputed from the input sequence
                intensities: new[] { 0.9, 0.4 });

            var predictions = model.TestResponseToPredictions(new[] { response }, inputs);

            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].FragmentAnnotations, Is.EquivalentTo(new[] { "b2+1", "y3+1" }));
            Assert.That(predictions[0].FragmentMZs!.All(mz => mz > 0), Is.True);
        }

        [Test]
        public void ResponseToPredictions_MapToInputFullSequence_DropsUnmappableAnnotations()
        {
            var model = new CoverageModel(FragmentIonMappingMode.MapToInputFullSequence);
            var inputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 30, null, null) { ValidatedFullSequence = "PEPTIDEK" }
            };
            // "b99" has no theoretical product for an 8-residue peptide and is dropped.
            var response = BuildJsonResponse(
                annotations: new[] { "b2+1", "b99+1" },
                mz: new[] { 0.0, 0.0 },
                intensities: new[] { 0.9, 0.4 });

            var predictions = model.TestResponseToPredictions(new[] { response }, inputs);

            Assert.That(predictions[0].FragmentAnnotations, Is.EquivalentTo(new[] { "b2+1" }));
        }

        [Test]
        public void ResponseToPredictions_MapToInputFullSequence_HandlesCaretChargeAnnotations()
        {
            // Altimeter annotations use '^' for charge; the old hardcoded '+' parsing threw here.
            var model = new TestableAltimeter(FragmentIonMappingMode.MapToInputFullSequence);
            var inputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 30, null, null) { ValidatedFullSequence = "PEPTIDEK" }
            };
            var response = BuildJsonResponse(
                annotations: new[] { "y3^1", "b2^1" },
                mz: new[] { 0.0, 0.0 },
                intensities: new[] { 0.7, 0.5 });

            var predictions = model.TestResponseToPredictions(new[] { response }, inputs);

            Assert.That(predictions[0].FragmentAnnotations, Is.EquivalentTo(new[] { "y3^1", "b2^1" }));
            Assert.That(predictions[0].FragmentMZs!.All(mz => mz > 0), Is.True);
        }

        [Test]
        public void ResponseToPredictions_MapToInputFullSequence_SubtractsNeutralLossMass()
        {
            var model = new TestableAltimeter(FragmentIonMappingMode.MapToInputFullSequence);
            var inputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 30, null, null) { ValidatedFullSequence = "PEPTIDEK" }
            };
            var response = BuildJsonResponse(
                annotations: new[] { "y3^1", "y3-H2O^1" },
                mz: new[] { 0.0, 0.0 },
                intensities: new[] { 0.7, 0.5 });

            var predictions = model.TestResponseToPredictions(new[] { response }, inputs);

            Assert.That(predictions[0].FragmentMZs!.Count, Is.EqualTo(2));
            var noLoss = predictions[0].FragmentMZs![0];
            var withLoss = predictions[0].FragmentMZs![1];
            Assert.That(noLoss - withLoss,
                Is.EqualTo(ChemicalFormula.ParseFormula("H2O").MonoisotopicMass).Within(0.001));
        }

        // ── GenerateLibrarySpectraFromPredictions ───────────────────────────────────

        [Test]
        public void GenerateLibrarySpectra_BuildsSpectrumWithMatchedPeaks()
        {
            var model = new CoverageModel();
            var prediction = new PeptideFragmentIntensityPrediction(
                "PEPTIDEK", "PEPTIDEK", 2,
                FragmentAnnotations: new List<string> { "b2+1", "y3+1" },
                FragmentMZs: new List<double> { 227.0, 375.0 },
                FragmentIntensities: new List<double> { 0.9, 0.4 });
            model.Seed(new List<PeptideFragmentIntensityPrediction> { prediction }, new[] { true });

            var spectra = model.GenerateLibrarySpectraFromPredictions(new double?[] { 30.0 }, out var warning);

            Assert.That(spectra.Count, Is.EqualTo(1));
            Assert.That(spectra[0].MatchedFragmentIons.Count, Is.EqualTo(2));
            Assert.That(spectra[0].ChargeState, Is.EqualTo(2));
            // No output filepath was provided, so the spectra-not-saved warning is expected.
            Assert.That(warning, Is.Not.Null);
            Assert.That(warning!.Message, Does.Contain("No file path"));
        }

        [Test]
        public void GenerateLibrarySpectra_FiltersLowIntensityAndImpossibleIons()
        {
            var model = new CoverageModel();
            var prediction = new PeptideFragmentIntensityPrediction(
                "PEPTIDEK", "PEPTIDEK", 2,
                FragmentAnnotations: new List<string> { "b2+1", "y3+1", "y4+1" },
                FragmentMZs: new List<double> { 227.0, 375.0, 400.0 },
                FragmentIntensities: new List<double> { 0.9, -1.0, 1e-9 });
            model.Seed(new List<PeptideFragmentIntensityPrediction> { prediction }, new[] { true });

            var spectra = model.GenerateLibrarySpectraFromPredictions(new double?[] { 30.0 }, out _, minIntensityFilter: 1e-4);

            // -1 (impossible) and the 1e-9 peak (below filter) are excluded; only b2 survives.
            Assert.That(spectra[0].MatchedFragmentIons.Count, Is.EqualTo(1));
        }

        [Test]
        public void GenerateLibrarySpectra_DeduplicatesByName()
        {
            var model = new CoverageModel();
            var prediction = new PeptideFragmentIntensityPrediction(
                "PEPTIDEK", "PEPTIDEK", 2,
                FragmentAnnotations: new List<string> { "b2+1" },
                FragmentMZs: new List<double> { 227.0 },
                FragmentIntensities: new List<double> { 0.9 });
            model.Seed(
                new List<PeptideFragmentIntensityPrediction> { prediction, prediction },
                new[] { true, true });

            var spectra = model.GenerateLibrarySpectraFromPredictions(new double?[] { 30.0, 30.0 }, out var warning);

            Assert.That(spectra.Count, Is.EqualTo(1));
            Assert.That(warning, Is.Not.Null);
            Assert.That(warning!.Message, Does.Contain("Duplicate"));
        }

        [Test]
        public void GenerateLibrarySpectra_EmptyPredictions_ReturnsEmpty()
        {
            var model = new CoverageModel();
            model.Seed(new List<PeptideFragmentIntensityPrediction>(), Array.Empty<bool>());

            var spectra = model.GenerateLibrarySpectraFromPredictions(Array.Empty<double?>(), out _);

            Assert.That(spectra, Is.Empty);
        }

        [Test]
        public void GenerateLibrarySpectra_RetentionTimeCountMismatch_Throws()
        {
            var model = new CoverageModel();
            var prediction = new PeptideFragmentIntensityPrediction(
                "PEPTIDEK", "PEPTIDEK", 2,
                new List<string> { "b2+1" }, new List<double> { 227.0 }, new List<double> { 0.9 });
            model.Seed(new List<PeptideFragmentIntensityPrediction> { prediction }, new[] { true });

            Assert.Throws<ArgumentException>(
                () => model.GenerateLibrarySpectraFromPredictions(new double?[] { 1.0, 2.0 }, out _));
        }

        [Test]
        public void GenerateLibrarySpectra_MapToValidatedFullSequence_LabelsSpectrumWithMassSourceSequence()
        {
            // When mod handling rewrites an incompatible modification, the validated (cleaned) sequence
            // differs chemically from the requested FullSequence. Under MapToValidatedFullSequence the
            // precursor m/z and peaks are built from the validated sequence, so the library label must
            // follow it too (peptide.FullSequence), not the original FullSequence — otherwise the
            // library identity desynchronizes from its spectral content.
            var model = new CoverageModel(FragmentIonMappingMode.MapToValidatedFullSequence);
            var prediction = new PeptideFragmentIntensityPrediction(
                FullSequence: "PEPM[Common Variable:Oxidation on M]TIDEK", // requested (modified)
                ValidatedFullSequence: "PEPMTIDEK",                        // cleaned: incompatible mod stripped
                PrecursorCharge: 2,
                FragmentAnnotations: new List<string> { "b2+1", "y3+1" },
                FragmentMZs: new List<double> { 227.0, 375.0 },
                FragmentIntensities: new List<double> { 0.9, 0.4 });
            model.Seed(new List<PeptideFragmentIntensityPrediction> { prediction }, new[] { true });

            var spectra = model.GenerateLibrarySpectraFromPredictions(new double?[] { 30.0 }, out _);

            Assert.That(spectra.Count, Is.EqualTo(1));
            Assert.That(spectra[0].Sequence, Is.EqualTo("PEPMTIDEK"),
                "Spectrum label must match the sequence the masses were built from (validated), not the requested FullSequence.");
        }

        [Test]
        public void GenerateLibrarySpectra_NeutralLossFragment_CarriesLossOnProduct()
        {
            // A y3 and its water-loss y3-H2O must become two distinct peaks. The loss peak's
            // theoretical product has to carry the loss so its annotation shows it and its mass
            // lines up with the loss-shifted m/z; reusing the base y3 product left the loss peak
            // labeled "y3" with a mass error of a whole water.
            var model = new CoverageModel();
            var prediction = new PeptideFragmentIntensityPrediction(
                "PEPTIDEK", "PEPTIDEK", 2,
                FragmentAnnotations: new List<string> { "y3+1", "y3-H2O+1" },
                FragmentMZs: new List<double> { 0.0, 0.0 },   // recomputed from theoretical products
                FragmentIntensities: new List<double> { 0.7, 0.5 });
            model.Seed(new List<PeptideFragmentIntensityPrediction> { prediction }, new[] { true });

            var spectra = model.GenerateLibrarySpectraFromPredictions(new double?[] { 30.0 }, out _);

            var peaks = spectra[0].MatchedFragmentIons;
            Assert.That(peaks.Count, Is.EqualTo(2));

            var waterMass = ChemicalFormula.ParseFormula("H2O").MonoisotopicMass;
            var lossPeak = peaks.Single(p => p.NeutralTheoreticalProduct.NeutralLoss != 0);
            var basePeak = peaks.Single(p => p.NeutralTheoreticalProduct.NeutralLoss == 0);

            Assert.That(lossPeak.NeutralTheoreticalProduct.NeutralLoss, Is.EqualTo(waterMass).Within(0.001));
            Assert.That(lossPeak.Annotation, Is.Not.EqualTo(basePeak.Annotation));
            // Product mass and m/z agree (no whole-water mass error), and the loss is the water shift.
            Assert.That(Math.Abs(lossPeak.MassErrorDa), Is.LessThan(0.001));
            Assert.That(basePeak.Mz - lossPeak.Mz, Is.EqualTo(waterMass).Within(0.001));
        }

        private sealed class CoverageModel : FragmentIntensityModel
        {
            private static readonly ISequenceConverter Conv = CreateUnimodConverter(
                UnimodSequenceFormatSchema.Instance, new HashSet<int> { 35, 4 });

            public override string ModelName => "CoverageModel";
            public override int MaxBatchSize => 1000;
            public override int MaxNumberOfBatchesPerRequest { get; init; } = 1;
            public override int ThrottlingDelayInMilliseconds { get; init; } = 0;
            public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 0;
            public override int MaxPeptideLength => 30;
            public override int MinPeptideLength => 1;
            public override HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6 };
            public override SequenceConversionHandlingMode ModHandlingMode { get; init; }
            public override IncompatibleParameterHandlingMode ParameterHandlingMode { get; init; }
            public override FragmentIonMappingMode FragmentIonMappingMode { get; init; }

            public CoverageModel(FragmentIonMappingMode mappingMode = FragmentIonMappingMode.MapToValidatedFullSequence)
                : base(Conv)
            {
                FragmentIonMappingMode = mappingMode;
                ModHandlingMode = SequenceConversionHandlingMode.ReturnNull;
                ParameterHandlingMode = IncompatibleParameterHandlingMode.ReturnNull;
            }

            protected override List<Dictionary<string, object>> ToBatchedRequests(List<FragmentIntensityPredictionInput> validInputs)
                => new();

            public void Seed(List<PeptideFragmentIntensityPrediction> predictions, bool[] mask)
            {
                Predictions = predictions;
                ValidInputsMask = mask;
            }

            public List<PeptideFragmentIntensityPrediction> TestResponseToPredictions(
                IReadOnlyList<string> responses, List<FragmentIntensityPredictionInput> requestInputs)
            {
                ModelInputs = requestInputs;
                return ResponseToPredictions(responses, requestInputs);
            }
        }

        // Real Altimeter parser (caret-charge annotations) exposed for the MapToInputFullSequence path.
        private sealed class TestableAltimeter : Altimeter2024Intensities
        {
            public TestableAltimeter(FragmentIonMappingMode mappingMode)
                : base(fragmentIonMappingMode: mappingMode) { }

            public List<PeptideFragmentIntensityPrediction> TestResponseToPredictions(
                IReadOnlyList<string> responses, List<FragmentIntensityPredictionInput> requestInputs)
            {
                ModelInputs = requestInputs;
                return ResponseToPredictions(responses, requestInputs);
            }
        }
    }
}
