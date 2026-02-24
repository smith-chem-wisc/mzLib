using NUnit.Framework;
using Readers.SpectralLibrary;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.IO;
using System;
using System.ComponentModel;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.AbstractClasses;

namespace Test.KoinaTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class Prosit2020IntensityHCDTests
    {
        [Test]
        public static void TestKoinaProsit2020IntensityHCDModelWritesReadableSpectralLibrary()
        {
            var experPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"SpectralLibrary\SpectralLibraryData\myPrositLib.msp");
            var predPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"SpectralLibrary\SpectralLibraryData\koinaTestOutput.msp");
            SpectralLibrary testLibraryWithoutDecoy = null;
            SpectralLibrary spectralLibraryTest = null;
            try
            {
                testLibraryWithoutDecoy = new SpectralLibrary(new List<string> { experPath });
                string aminoacids = @"ACDEFGHIKLMNPQRSTVWY";
                var librarySpectra = testLibraryWithoutDecoy.GetAllLibrarySpectra()
                    .Where(p => p.ChargeState < 6 && p.Sequence.Length < 30 && p.Sequence.Length > 1 && p.Sequence.ToHashSet().IsSubsetOf(aminoacids.ToHashSet())).ToList();

                var peptides = librarySpectra.Select(p => p.Sequence).ToList();
                var charges = librarySpectra.Select(p => p.ChargeState).ToList();
                var energies = librarySpectra.Select(p => 35).ToList();
                var retentionTimes = librarySpectra.Select(p => p.RetentionTime).ToArray();

                // Create model inputs
                var modelInputs = new List<FragmentIntensityPredictionInput>();
                for (int i = 0; i < peptides.Count; i++)
                {
                    modelInputs.Add(new FragmentIntensityPredictionInput(
                        FullSequence: peptides[i],
                        PrecursorCharge: charges[i],
                        CollisionEnergy: energies[i],
                        InstrumentType: null,
                        FragmentationType: null
                    ));
                }

                var modelHandler = new Prosit2020IntensityHCD();
                var predictions = modelHandler.Predict(modelInputs);
                
                // Generate spectral library from predictions
                var predictedSpectra = modelHandler.GenerateLibrarySpectraFromPredictions(retentionTimes, out var warning, predPath, minIntensityFilter: 1e-6);

                // Test that the predicted spectrum that was saved matches the predicted spectra in memory
                spectralLibraryTest = new SpectralLibrary(new List<string> { predPath });
                var spectralLibrary = spectralLibraryTest.GetAllLibrarySpectra().ToList();
                Assert.That(peptides.Count == spectralLibrary.Count);

                var sortedPredictedSpectra = predictedSpectra.OrderBy(p => p.Sequence).ThenBy(p => p.ChargeState).ToList();
                var sortedSavedSpectra = spectralLibrary.OrderBy(p => p.Sequence).ThenBy(p => p.ChargeState).ToList();
                for (int i = 0; i < peptides.Count; i++)
                {
                    var inMemorySpectrum = sortedPredictedSpectra[i];
                    var savedSpectrum = sortedSavedSpectra[i];
                    Assert.That(inMemorySpectrum.Sequence == savedSpectrum.Sequence);
                    Assert.That(inMemorySpectrum.ChargeState == savedSpectrum.ChargeState);
                    Assert.That(inMemorySpectrum.MatchedFragmentIons.Count == savedSpectrum.MatchedFragmentIons.Count);

                    var sortedInMemoryFrags = inMemorySpectrum.MatchedFragmentIons.OrderBy(p => p.Mz).ToList();
                    var sortedSavedFrags = savedSpectrum.MatchedFragmentIons.OrderBy(p => p.Mz).ToList();

                    // Get max intensity for normalization
                    double maxInMemoryIntensity = sortedInMemoryFrags.Max(f => f.Intensity);

                    for (int j = 0; j < inMemorySpectrum.MatchedFragmentIons.Count; j++)
                    {
                        var inMemoryFrag = sortedInMemoryFrags[j];
                        var savedFrag = sortedSavedFrags[j];
                        Assert.That(inMemoryFrag.Mz, Is.EqualTo(savedFrag.Mz).Within(1e-6));
                        // Compare normalized intensities since saved library normalizes to max intensity
                        Assert.That(inMemoryFrag.Intensity / maxInMemoryIntensity, Is.EqualTo(savedFrag.Intensity).Within(1e-6));
                        Assert.That(inMemoryFrag.NeutralTheoreticalProduct.ProductType == savedFrag.NeutralTheoreticalProduct.ProductType);
                        Assert.That(inMemoryFrag.NeutralTheoreticalProduct.FragmentNumber == savedFrag.NeutralTheoreticalProduct.FragmentNumber);
                        Assert.That(inMemoryFrag.Charge == savedFrag.Charge);
                    }
                }
            }
            finally
            {
                testLibraryWithoutDecoy?.CloseConnections();
                spectralLibraryTest?.CloseConnections();
                File.Delete(predPath);
            }
        }

        /// <summary>
        /// Tests that the Prosit2020IntensityHCD model accepts valid modified peptide sequences.
        /// </summary>
        [Test]
        public static void TestKoinaProsit2020IntensityHCDModelAcceptsValidPeptidesWithModifications()
        {
            // Valid peptide with modifications (string length > 30 with mods, but valid base sequence length)
            var modelInputs = new List<FragmentIntensityPredictionInput>
            {
                new FragmentIntensityPredictionInput(
                    FullSequence: "M[Common Variable:Oxidation on M]PEPTIDEC[Common Fixed:Carbamidomethyl on C]A",
                    PrecursorCharge: 2,
                    CollisionEnergy: 35,
                    InstrumentType: null,
                    FragmentationType: null
                ),
                new FragmentIntensityPredictionInput(
                    FullSequence: "LMNPQRSTVWY",
                    PrecursorCharge: 2,
                    CollisionEnergy: 35,
                    InstrumentType: null,
                    FragmentationType: null
                )
            };

            var model = new Prosit2020IntensityHCD();
            var predictions = model.Predict(modelInputs);
            
            Assert.That(predictions.Count, Is.EqualTo(2));
            Assert.That(predictions.All(p => p.Warning == null || p.Warning.Message.Contains("skipped")));
        }

        /// <summary>
        /// Tests that empty input lists are handled gracefully without throwing exceptions.
        /// This is a common edge case when filtering results in no valid peptides.
        /// </summary>
        [Test]
        public static void TestKoinaProsit2020IntensityHCDModelEmptyInputHandling()
        {
            var emptyInputs = new List<FragmentIntensityPredictionInput>();

            var model = new Prosit2020IntensityHCD();
            var predictions = model.Predict(emptyInputs);

            Assert.That(predictions.Count, Is.EqualTo(0));
            Assert.DoesNotThrow( () => model.Predict(emptyInputs));
        }

        /// <summary>
        /// Tests boundary conditions for charge states.
        /// Valid charge states are 1-6 inclusive for the Prosit 2020 model.
        /// </summary>
        [Test]
        public static void TestKoinaProsit2020IntensityHCDModelChargeStateBoundaries()
        {
            var peptide = "PEPTIDEK";
            var energy = 35;

            var model = new Prosit2020IntensityHCD();

            // Test all valid charge states (1-6)
            for (int charge = 1; charge <= 6; charge++)
            {
                var modelInputs = new List<FragmentIntensityPredictionInput>
                {
                    new FragmentIntensityPredictionInput(
                        FullSequence: peptide,
                        PrecursorCharge: charge,
                        CollisionEnergy: energy,
                        InstrumentType: null,
                        FragmentationType: null
                    )
                };

                var predictions = model.Predict(modelInputs);

                Assert.That(predictions.Count, Is.EqualTo(1),
                    $"Charge state {charge} should be valid");
                Assert.That(predictions[0].Warning, Is.Null,
                    $"Charge state {charge} should not produce a warning");
            }

            // Test invalid charge state 0
            var modelInputsZero = new List<FragmentIntensityPredictionInput>
            {
                new FragmentIntensityPredictionInput(
                    FullSequence: peptide,
                    PrecursorCharge: 0,
                    CollisionEnergy: energy,
                    InstrumentType: null,
                    FragmentationType: null
                )
            };

            var predictionsZero = model.Predict(modelInputsZero);

            Assert.That(predictionsZero.Count, Is.EqualTo(1),
                "Should return a prediction entry for invalid charge");
            Assert.That(predictionsZero[0].Warning, Is.Not.Null,
                "Charge state 0 should produce a warning");
            Assert.That(predictionsZero[0].FragmentAnnotations, Is.Null,
                "Invalid charge should result in null predictions");

            // Test invalid negative charge state
            var modelInputsNegative = new List<FragmentIntensityPredictionInput>
            {
                new FragmentIntensityPredictionInput(
                    FullSequence: peptide,
                    PrecursorCharge: -1,
                    CollisionEnergy: energy,
                    InstrumentType: null,
                    FragmentationType: null
                )
            };

            var predictionsNegative = model.Predict(modelInputsNegative);

            Assert.That(predictionsNegative.Count, Is.EqualTo(1),
                "Should return a prediction entry for negative charge");
            Assert.That(predictionsNegative[0].Warning, Is.Not.Null,
                "Negative charge state should produce a warning");
        }

        /// <summary>
        /// Tests boundary conditions for peptide sequence length.
        /// Valid lengths are 1-30 characters (after removing modification annotations).
        /// </summary>
        [Test]
        public static void TestKoinaProsit2020IntensityHCDModelPeptideLengthBoundaries()
        {
            var charge = 2;
            var energy = 35;
            var model = new Prosit2020IntensityHCD();

            // Test minimum valid length (1 amino acid)
            var minLengthInputs = new List<FragmentIntensityPredictionInput>
            {
                new FragmentIntensityPredictionInput(
                    FullSequence: "K",
                    PrecursorCharge: charge,
                    CollisionEnergy: energy,
                    InstrumentType: null,
                    FragmentationType: null
                )
            };

            var predictionsMin = model.Predict(minLengthInputs);

            Assert.That(predictionsMin.Count, Is.EqualTo(1),
                "Single amino acid should be valid");
            Assert.That(predictionsMin[0].Warning, Is.Null,
                "Minimum length peptide should not produce a warning");

            // Test maximum valid length (30 amino acids)
            var maxLengthInputs = new List<FragmentIntensityPredictionInput>
            {
                new FragmentIntensityPredictionInput(
                    FullSequence: new string('A', 30),
                    PrecursorCharge: charge,
                    CollisionEnergy: energy,
                    InstrumentType: null,
                    FragmentationType: null
                )
            };

            var predictionsMax = model.Predict(maxLengthInputs);

            Assert.That(predictionsMax.Count, Is.EqualTo(1),
                "Peptide with exactly 30 amino acids should be valid");
            Assert.That(predictionsMax[0].Warning, Is.Null,
                "Maximum length peptide should not produce a warning");

            // Test invalid length (0 amino acids)
            var tooShortInputs = new List<FragmentIntensityPredictionInput>
            {
                new FragmentIntensityPredictionInput(
                    FullSequence: "",
                    PrecursorCharge: charge,
                    CollisionEnergy: energy,
                    InstrumentType: null,
                    FragmentationType: null
                )
            };

            var predictionsTooShort = model.Predict(tooShortInputs);

            Assert.That(predictionsTooShort.Count, Is.EqualTo(1),
                "Should return a prediction entry for empty sequence");
            Assert.That(predictionsTooShort[0].Warning, Is.Not.Null,
                "Too short peptide should produce a warning");

            // Test invalid length (31 amino acids)
            var tooLongInputs = new List<FragmentIntensityPredictionInput>
            {
                new FragmentIntensityPredictionInput(
                    FullSequence: new string('A', 31),
                    PrecursorCharge: charge,
                    CollisionEnergy: energy,
                    InstrumentType: null,
                    FragmentationType: null
                )
            };

            var predictionsTooLong = model.Predict(tooLongInputs);

            Assert.That(predictionsTooLong.Count, Is.EqualTo(1),
                "Should return a prediction entry for too long sequence");
            Assert.That(predictionsTooLong[0].Warning, Is.Not.Null,
                "Too long peptide should produce a warning");
        }

        /// <summary>
        /// Tests handling of peptides containing invalid or non-standard amino acids.
        /// The standard 20 amino acids are: ACDEFGHIKLMNPQRSTVWY
        /// Characters like X (unknown), B (Asx), Z (Glx), U (Selenocysteine), O (Pyrrolysine)
        /// may not be supported by the model.
        /// </summary>
        [Test]
        public static void TestKoinaProsit2020IntensityHCDModelInvalidAminoAcids()
        {
            var charge = 2;
            var energy = 35;
            var model = new Prosit2020IntensityHCD();

            // Test peptide with 'X' (unknown amino acid) - common in database searches
            var inputsWithX = new List<FragmentIntensityPredictionInput>
            {
                new FragmentIntensityPredictionInput(
                    FullSequence: "PEPTXIDEK",
                    PrecursorCharge: charge,
                    CollisionEnergy: energy,
                    InstrumentType: null,
                    FragmentationType: null
                )
            };

            var predictionsX = model.Predict(inputsWithX);

            // Document expected behavior - should either filter out or handle gracefully
            Assert.That(predictionsX.Count, Is.EqualTo(1),
                "Should return a prediction entry");
            Assert.That(predictionsX[0].FragmentAnnotations, Is.Null,
                "Peptide with 'X' (unknown AA) should be filtered out");

            // Test peptide with 'U' (Selenocysteine) - rare but valid amino acid
            var inputsWithU = new List<FragmentIntensityPredictionInput>
            {
                new FragmentIntensityPredictionInput(
                    FullSequence: "PEPTUIDEK",
                    PrecursorCharge: charge,
                    CollisionEnergy: energy,
                    InstrumentType: null,
                    FragmentationType: null
                )
            };

            var predictionsU = model.Predict(inputsWithU);

            // Selenocysteine is typically not supported by most models
            Assert.That(predictionsU.Count, Is.EqualTo(1),
                "Should return a prediction entry");
            Assert.That(predictionsU[0].FragmentAnnotations, Is.Null,
                "Peptide with 'U' (Selenocysteine) should be filtered out if not supported");

            // Test mixed valid and invalid - ensure valid ones are kept
            var mixedInputs = new List<FragmentIntensityPredictionInput>
            {
                new FragmentIntensityPredictionInput("PEPTIDEK", charge, energy, null, null),
                new FragmentIntensityPredictionInput("PEPTXIDEK", charge, energy, null, null),
                new FragmentIntensityPredictionInput("VALIDSEQ", charge, energy, null, null)
            };

            var predictionsMixed = model.Predict(mixedInputs);

            Assert.That(predictionsMixed.Count, Is.EqualTo(3),
                "Should return entries for all inputs");
            Assert.That(predictionsMixed.Count(p => p.FragmentAnnotations != null), Is.EqualTo(2),
                "Only valid peptides should have predictions");
        }

        /// <summary>
        /// Tests handling of duplicate peptide+charge combinations.
        /// Documents whether duplicates are deduplicated, kept, or cause errors.
        /// </summary>
        [Test]
        public static void TestKoinaProsit2020IntensityHCDModelDuplicateHandling()
        {
            var energy = 35;
            var model = new Prosit2020IntensityHCD();

            // Same peptide, same charge - exact duplicate
            var duplicateInputs = new List<FragmentIntensityPredictionInput>
            {
                new FragmentIntensityPredictionInput("PEPTIDEK", 2, energy, null, null),
                new FragmentIntensityPredictionInput("PEPTIDEK", 2, energy, null, null)
            };

            var predictionsDuplicate = model.Predict(duplicateInputs);

            // Document the expected behavior for duplicates
            Assert.That(predictionsDuplicate.Count, Is.EqualTo(2),
                "Duplicate peptides should be handled without error");

            // Same peptide, different charge - should both be kept (not duplicates)
            var differentChargeInputs = new List<FragmentIntensityPredictionInput>
            {
                new FragmentIntensityPredictionInput("PEPTIDEK", 2, energy, null, null),
                new FragmentIntensityPredictionInput("PEPTIDEK", 3, energy, null, null)
            };

            var predictionsDiffCharge = model.Predict(differentChargeInputs);

            Assert.That(predictionsDiffCharge.Count, Is.EqualTo(2),
                "Same peptide with different charges should both be kept");
            Assert.That(predictionsDiffCharge.All(p => p.FragmentAnnotations != null), Is.True,
                "Both predictions should be valid");
        }

        /// <summary>
        /// Tests that null retention times are handled correctly during spectral library generation.
        /// Retention times are optional metadata and may not always be available.
        /// </summary>
        [Test]
        public static void TestKoinaProsit2020IntensityHCDModelNullRetentionTimes()
        {
            var modelInputs = new List<FragmentIntensityPredictionInput>
            {
                new FragmentIntensityPredictionInput("PEPTIDEK", 2, 35, null, null),
                new FragmentIntensityPredictionInput("VALIDSEQK", 2, 35, null, null),
                new FragmentIntensityPredictionInput("ANTHERSEQR", 3, 35, null, null)
            };

            var model = new Prosit2020IntensityHCD();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(3),
                "All peptides should have predictions");

            // Test spectral library generation with mixed null/non-null retention times
            var retentionTimes = new double[] { 100.0, 0.0, 200.0 };
            Assert.DoesNotThrow(() => 
                model.GenerateLibrarySpectraFromPredictions(retentionTimes, out var warning),
                "Mixed retention times should be handled gracefully");
        }

        /// <summary>
        /// Tests that predicted spectra contain chemically reasonable values.
        /// Fragment m/z values should be positive and within expected mass range.
        /// Intensities should be non-negative.
        /// </summary>
        [Test]
        public static void TestKoinaProsit2020IntensityHCDModelPredictionQuality()
        {
            var modelInputs = new List<FragmentIntensityPredictionInput>
            {
                new FragmentIntensityPredictionInput("PEPTIDEK", 2, 35, null, null),
                new FragmentIntensityPredictionInput("ELVISLIVESK", 2, 35, null, null)
            };

            var model = new Prosit2020IntensityHCD();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(2),
                "Should have predictions for both peptides");

            foreach (var prediction in predictions.Where(p => p.FragmentAnnotations != null))
            {
                Assert.That(prediction.FragmentAnnotations.Count, Is.GreaterThan(0),
                    $"Spectrum for {prediction.FullSequence} should have fragment ions");

                for (int i = 0; i < prediction.FragmentAnnotations.Count; i++)
                {
                    // m/z should be positive and reasonable for peptide fragments
                    Assert.That(prediction.FragmentMZs[i], Is.GreaterThan(0),
                        "Fragment m/z should be positive");
                    Assert.That(prediction.FragmentMZs[i], Is.LessThan(5000),
                        "Fragment m/z should be within reasonable range for peptide fragments");

                    // Intensity should be non-negative (or -1 for impossible ions)
                    Assert.That(prediction.FragmentIntensities[i], Is.GreaterThanOrEqualTo(-1),
                        "Fragment intensity should be non-negative or -1 for impossible ions");
                }
            }
        }

        [Test]
        public static void TestKoinaProsit2020IntensityHCDModelFiltersDuplicatePeptidesFromSpectralLibrary()
        {
            var modelInputs = new List<FragmentIntensityPredictionInput>
            {
                new FragmentIntensityPredictionInput("PEPTIDEK", 2, 35, null, null),
                new FragmentIntensityPredictionInput("PEPTIDEK", 2, 35, null, null),
                new FragmentIntensityPredictionInput("ELVISLIVESK", 2, 35, null, null)
            };

            var model = new Prosit2020IntensityHCD();
            var predictions =  model.Predict(modelInputs);
            
            Assert.That(predictions.Count, Is.EqualTo(3),
                "All inputs should have prediction entries");

            var retentionTimes = new double[] { 100.0, 100.0, 200.0 };
            var spectra = model.GenerateLibrarySpectraFromPredictions(retentionTimes, out var warning);
            
            Assert.That(spectra.Count, Is.EqualTo(2),
                "Duplicate peptide should be filtered out during library generation");
            Assert.That(warning, Is.Not.Null,
                "Duplicate spectra warning should be generated when duplicates are present");
        }

        [Test]
        public static void TestKoinaProsit2020IntensityHCDModelRequestBatching()
        {
            var aminoacids = "ACDEFGHIKLMNPQRSTVWY".ToArray();
            var modelInputs = new List<FragmentIntensityPredictionInput>();
            var seqLength = 20;
            var numberOfSequences = 2500;

            var peptides = new HashSet<string>();
            while (peptides.Count < numberOfSequences)
            {
                var pep = new string(Random.Shared.GetItems(aminoacids, seqLength));
                peptides.Add(pep);
            }

            foreach (var peptide in peptides)
            {
                modelInputs.Add(new FragmentIntensityPredictionInput(peptide, 2, 35, null, null));
            }

            var model = new Prosit2020IntensityHCD();
            Assert.DoesNotThrow( () => model.Predict(modelInputs));
            Assert.That(model.Predictions.Count, Is.EqualTo(numberOfSequences));
        }
    }
}
