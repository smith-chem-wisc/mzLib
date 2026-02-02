using NUnit.Framework;
using Readers.SpectralLibrary;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.IO;
using System;
using System.ComponentModel;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;

namespace Test.KoinaTests.FragmentIntensityPrediction
{
    public class Koina_Prosit2020IntensityHCD
    {
        [Test]
        public static async Task TestKoinaProsit2020IntensityHCDModelWritesReadableSpectralLibrary()
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
                var retentionTimes = librarySpectra.Select(p => p.RetentionTime).ToList();

                var modelHandler = new Prosit2020IntensityHCD(peptides, charges, energies, retentionTimes, out var warnings, minIntensityFilter: 1e-6);

                await modelHandler.RunInferenceAsync();
                var predictedSpectra = modelHandler.PredictedSpectra;
                modelHandler.SavePredictedSpectralLibrary(predPath);

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
            catch
            {
                Assert.Fail("Test failed somewhere in the try statement.");
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
        public static void TestKoinaProsit2020IntensityHCDModelAcceptsValidPeptidesWithModifications()
        { 
            var validCharges = new List<int> { 2, 2 };
            var validEnergies = new List<int> { 35, 35 };
            var validRetentionTimes = new List<double?> { 100.0, 200.0 };
            WarningException warning;

            // Valid peptide with modifications (string length > 30 with mods, but valid base sequence length)
            var peptidesWithMods = new List<string>
            {
                "M[Common Variable:Oxidation on M]PEPTIDEC[Common Fixed:Carbamidomethyl on C]A",
                "LMNPQRSTVWY"
            };
            var model = new Prosit2020IntensityHCD(peptidesWithMods, validCharges, validEnergies, validRetentionTimes, out warning);
            Assert.That(warning, Is.Null);
            Assert.That(model.PeptideSequences.Count == 2);
        }

        /// <summary>
        /// Tests that the input lists need to be mappable. If lengths do not match, an exception is thrown.
        /// </summary>
        [Test]
        public static void TestKoinaProsit2020IntensityHCDModelThrowsOnMismatchedInputLengths()
        {
            var validPeptides = new List<string> { "ACDEFGHIK", "LMNPQRSTVWY" };
            var validEnergies = new List<int> { 35, 35 };
            var validRetentionTimes = new List<double?> { 100.0, 200.0 };

            // Mismatched lengths
            var mismatchedCharges = new List<int> { 2 };
            Assert.Throws<ArgumentException>(() =>
                new Prosit2020IntensityHCD(validPeptides, mismatchedCharges, validEnergies, validRetentionTimes, out var _));
        }

        /// <summary>
        /// Tests that empty input lists are handled gracefully without throwing exceptions.
        /// This is a common edge case when filtering results in no valid peptides.
        /// </summary>
        [Test]
        public static void TestKoinaProsit2020IntensityHCDModelEmptyInputHandling()
        {
            var emptyPeptides = new List<string>();
            var emptyCharges = new List<int>();
            var emptyEnergies = new List<int>();
            var emptyRetentionTimes = new List<double?>();

            // Empty lists should not throw - this can happen when all peptides are filtered out
            var model = new Prosit2020IntensityHCD(emptyPeptides, emptyCharges, emptyEnergies, emptyRetentionTimes, out var warning);

            Assert.That(model.PeptideSequences.Count, Is.EqualTo(0));
            Assert.That(warning, Is.Not.Null);
            Assert.DoesNotThrowAsync(async () => await model.RunInferenceAsync());
            Assert.That(model.PredictedSpectra.Count == 0);
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
            double? retentionTime = 100.0;

            // Test all valid charge states (1-6)
            for (int charge = 1; charge <= 6; charge++)
            {
                var model = new Prosit2020IntensityHCD(
                    new List<string> { peptide },
                    new List<int> { charge },
                    new List<int> { energy },
                    new List<double?> { retentionTime },
                    out var warning);

                Assert.That(model.PeptideSequences.Count, Is.EqualTo(1),
                    $"Charge state {charge} should be valid");
                Assert.That(warning, Is.Null,
                    $"Charge state {charge} should not produce a warning");
            }

            // Test invalid charge state 0
            var modelZero = new Prosit2020IntensityHCD(
                new List<string> { peptide },
                new List<int> { 0 },
                new List<int> { energy },
                new List<double?> { retentionTime },
                out var warningZero);

            Assert.That(modelZero.PeptideSequences.Count, Is.EqualTo(0),
                "Charge state 0 should be invalid");
            Assert.That(warningZero, Is.Not.Null,
                "Charge state 0 should produce a warning");

            // Test invalid negative charge state
            var modelNegative = new Prosit2020IntensityHCD(
                new List<string> { peptide },
                new List<int> { -1 },
                new List<int> { energy },
                new List<double?> { retentionTime },
                out var warningNegative);

            Assert.That(modelNegative.PeptideSequences.Count, Is.EqualTo(0),
                "Negative charge state should be invalid");
            Assert.That(warningNegative, Is.Not.Null,
                "Negative charge state should produce a warning");
        }

        /// <summary>
        /// Tests boundary conditions for peptide sequence length.
        /// Valid lengths are 2-30 characters (after removing modification annotations).
        /// </summary>
        [Test]
        public static void TestKoinaProsit2020IntensityHCDModelPeptideLengthBoundaries()
        {
            var charge = 2;
            var energy = 35;
            double? retentionTime = 100.0;

            // Test minimum valid length (2 amino acids)
            var minLengthPeptide = "KR";
            var modelMin = new Prosit2020IntensityHCD(
                new List<string> { minLengthPeptide },
                new List<int> { charge },
                new List<int> { energy },
                new List<double?> { retentionTime },
                out var warningMin);

            Assert.That(modelMin.PeptideSequences.Count, Is.EqualTo(1),
                "Peptide with 2 amino acids should be valid");
            Assert.That(warningMin, Is.Null,
                "Minimum length peptide should not produce a warning");

            // Test maximum valid length (30 amino acids)
            var maxLengthPeptide = new string('A', 30);
            var modelMax = new Prosit2020IntensityHCD(
                new List<string> { maxLengthPeptide },
                new List<int> { charge },
                new List<int> { energy },
                new List<double?> { retentionTime },
                out var warningMax);

            Assert.That(modelMax.PeptideSequences.Count, Is.EqualTo(1),
                "Peptide with exactly 30 amino acids should be valid");
            Assert.That(warningMax, Is.Null,
                "Maximum length peptide should not produce a warning");

            // Test invalid length (0 amino acid)
            var tooShortPeptide = "";
            var modelTooShort = new Prosit2020IntensityHCD(
                new List<string> { tooShortPeptide },
                new List<int> { charge },
                new List<int> { energy },
                new List<double?> { retentionTime },
                out var warningTooShort);

            Assert.That(modelTooShort.PeptideSequences.Count, Is.EqualTo(0),
                "Peptide with 0 amino acids should be invalid");
            Assert.That(warningTooShort, Is.Not.Null,
                "Too short peptide should produce a warning");

            // Test invalid length (31 amino acids)
            var tooLongPeptide = new string('A', 31);
            var modelTooLong = new Prosit2020IntensityHCD(
                new List<string> { tooLongPeptide },
                new List<int> { charge },
                new List<int> { energy },
                new List<double?> { retentionTime },
                out var warningTooLong);

            Assert.That(modelTooLong.PeptideSequences.Count, Is.EqualTo(0),
                "Peptide with 31 amino acids should be invalid");
            Assert.That(warningTooLong, Is.Not.Null,
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
            double? retentionTime = 100.0;

            // Test peptide with 'X' (unknown amino acid) - common in database searches
            var peptideWithX = "PEPTXIDEK";
            var modelX = new Prosit2020IntensityHCD(
                new List<string> { peptideWithX },
                new List<int> { charge },
                new List<int> { energy },
                new List<double?> { retentionTime },
                out var warningX);

            // Document expected behavior - should either filter out or handle gracefully
            Assert.That(modelX.PeptideSequences.Count, Is.EqualTo(0),
                "Peptide with 'X' (unknown AA) should be filtered out");

            // Test peptide with 'U' (Selenocysteine) - rare but valid amino acid
            var peptideWithU = "PEPTUIDEK";
            var modelU = new Prosit2020IntensityHCD(
                new List<string> { peptideWithU },
                new List<int> { charge },
                new List<int> { energy },
                new List<double?> { retentionTime },
                out var warningU);

            // Selenocysteine is typically not supported by most models
            Assert.That(modelU.PeptideSequences.Count, Is.EqualTo(0),
                "Peptide with 'U' (Selenocysteine) should be filtered out if not supported");

            // Test mixed valid and invalid - ensure valid ones are kept
            var mixedPeptides = new List<string> { "PEPTIDEK", "PEPTXIDEK", "VALIDSEQ" };
            var mixedCharges = new List<int> { charge, charge, charge };
            var mixedEnergies = new List<int> { energy, energy, energy };
            var mixedRetentionTimes = new List<double?> { retentionTime, retentionTime, retentionTime };

            var modelMixed = new Prosit2020IntensityHCD(
                mixedPeptides, mixedCharges, mixedEnergies, mixedRetentionTimes, out var warningMixed);

            Assert.That(modelMixed.PeptideSequences.Count, Is.EqualTo(2),
                "Only valid peptides should be retained when mixed with invalid ones");
        }

        /// <summary>
        /// Tests handling of duplicate peptide+charge combinations.
        /// Documents whether duplicates are deduplicated, kept, or cause errors.
        /// </summary>
        [Test]
        public static void TestKoinaProsit2020IntensityHCDModelDuplicateHandling()
        {
            var energy = 35;
            double? retentionTime = 100.0;

            // Same peptide, same charge - exact duplicate
            var duplicatePeptides = new List<string> { "PEPTIDEK", "PEPTIDEK" };
            var duplicateCharges = new List<int> { 2, 2 };
            var duplicateEnergies = new List<int> { energy, energy };
            var duplicateRetentionTimes = new List<double?> { retentionTime, retentionTime };

            var modelDuplicate = new Prosit2020IntensityHCD(
                duplicatePeptides, duplicateCharges, duplicateEnergies, duplicateRetentionTimes,
                out var warningDuplicate);

            // Document the expected behavior for duplicates
            // Option 1: Duplicates are kept (count == 2)
            // Option 2: Duplicates are deduplicated (count == 1)
            Assert.That(modelDuplicate.PeptideSequences.Count, Is.GreaterThan(0),
                "Duplicate peptides should be handled without error");

            // Same peptide, different charge - should both be kept (not duplicates)
            var samePeptideDiffCharge = new List<string> { "PEPTIDEK", "PEPTIDEK" };
            var differentCharges = new List<int> { 2, 3 };

            var modelDiffCharge = new Prosit2020IntensityHCD(
                samePeptideDiffCharge, differentCharges, duplicateEnergies, duplicateRetentionTimes,
                out var warningDiffCharge);

            Assert.That(modelDiffCharge.PeptideSequences.Count, Is.EqualTo(2),
                "Same peptide with different charges should both be kept");
            Assert.That(warningDiffCharge, Is.Null,
                "Same peptide with different charges should not produce a warning");
        }

        /// <summary>
        /// Tests that null retention times are handled correctly.
        /// Retention times are optional metadata and may not always be available.
        /// </summary>
        [Test]
        public static void TestKoinaProsit2020IntensityHCDModelNullRetentionTimes()
        {
            var peptides = new List<string> { "PEPTIDEK", "VALIDSEQK", "ANTHERSEQR" };
            var charges = new List<int> { 2, 2, 3 };
            var energies = new List<int> { 35, 35, 35 };

            // Mix of null and non-null retention times
            var retentionTimes = new List<double?> { 100.0, null, 200.0 };

            var model = new Prosit2020IntensityHCD(
                peptides, charges, energies, retentionTimes, out var warning);

            Assert.That(model.PeptideSequences.Count, Is.EqualTo(3),
                "Null retention times should not cause peptides to be filtered out");
            Assert.That(warning, Is.Null,
                "Null retention times should not produce a warning");

            // All null retention times
            var allNullRetentionTimes = new List<double?> { null, null, null };

            var modelAllNull = new Prosit2020IntensityHCD(
                peptides, charges, energies, allNullRetentionTimes, out var warningAllNull);

            Assert.That(modelAllNull.PeptideSequences.Count, Is.EqualTo(3),
                "All null retention times should still allow predictions");
        }

        /// <summary>
        /// Tests that predicted spectra contain chemically reasonable values.
        /// Fragment m/z values should be positive and within expected mass range.
        /// Intensities should be non-negative.
        /// </summary>
        [Test]
        public static async Task TestKoinaProsit2020IntensityHCDModelPredictionQuality()
        {
            var peptides = new List<string> { "PEPTIDEK", "ELVISLIVESK" };
            var charges = new List<int> { 2, 2 };
            var energies = new List<int> { 35, 35 };
            var retentionTimes = new List<double?> { 100.0, 200.0 };

            var model = new Prosit2020IntensityHCD(peptides, charges, energies, retentionTimes, out var warning);
            await model.RunInferenceAsync();

            Assert.That(model.PredictedSpectra.Count, Is.EqualTo(2),
                "Should have predictions for both peptides");

            foreach (var spectrum in model.PredictedSpectra)
            {
                Assert.That(spectrum.MatchedFragmentIons.Count, Is.GreaterThan(0),
                    $"Spectrum for {spectrum.Sequence} should have fragment ions");

                foreach (var fragment in spectrum.MatchedFragmentIons)
                {
                    // m/z should be positive and reasonable for peptide fragments
                    Assert.That(fragment.Mz, Is.GreaterThan(0),
                        "Fragment m/z should be positive");
                    Assert.That(fragment.Mz, Is.LessThan(5000),
                        "Fragment m/z should be within reasonable range for peptide fragments");

                    // Intensity should be non-negative
                    Assert.That(fragment.Intensity, Is.GreaterThanOrEqualTo(0),
                        "Fragment intensity should be non-negative");

                    // Charge should be positive
                    Assert.That(fragment.Charge, Is.GreaterThan(0),
                        "Fragment charge should be positive");
                }
            }
        }

        [Test]
        public static void TestKoinaProsit2020IntensityHCDModelRequestBatching()
        {
            var aminoacids = "ACDEFGHIKLMNPQRSTVWY".ToArray();
            var peptides = new List<string>();
            var charges = new List<int>();
            var energies = new List<int>();
            var retentionTimes = new List<double?>();
            var seqLength = 20;
            var numberOfSequences = 2500;

            while (peptides.Count < numberOfSequences)
            {
                var pep = new string(Random.Shared.GetItems(aminoacids, seqLength));
                if (!peptides.Contains(pep))
                {
                    peptides.Add(pep);
                    charges.Add(2);
                    energies.Add(35);
                    retentionTimes.Add(0);
                }
            }
            var modelHandler = new Prosit2020IntensityHCD(peptides, charges, energies, retentionTimes, out var warnings);
            Assert.DoesNotThrowAsync(async () => await modelHandler.RunInferenceAsync());
            Assert.That(modelHandler.PredictedSpectra.Count == numberOfSequences);
        }
    }
}
