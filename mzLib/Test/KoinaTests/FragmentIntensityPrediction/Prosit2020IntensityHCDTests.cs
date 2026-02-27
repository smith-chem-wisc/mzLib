using NUnit.Framework;
using Readers.SpectralLibrary;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.IO;
using System;
using System.ComponentModel;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;

namespace Test.KoinaTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class Prosit2020IntensityHCDTests
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

                var modelHandler = new Prosit2020IntensityHCD(peptides, charges, energies, retentionTimes, out var warnings, predPath, minIntensityFilter: 1e-6);

                WarningException? duplicatesWarning = await modelHandler.RunInferenceAsync();
                var predictedSpectra = modelHandler.PredictedSpectra;

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

            // Test minimum valid length (1 amino acids)
            var minLengthPeptide = "K";
            var modelMin = new Prosit2020IntensityHCD(
                new List<string> { minLengthPeptide },
                new List<int> { charge },
                new List<int> { energy },
                new List<double?> { retentionTime },
                out var warningMin);

            Assert.That(modelMin.PeptideSequences.Count, Is.EqualTo(1),
                "Single amino acid should be valid");
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
            WarningException? duplicatesWarning = await model.RunInferenceAsync();

            Assert.That(model.PredictedSpectra.Count, Is.EqualTo(2),
                "Should have predictions for both peptides");
            Assert.That(duplicatesWarning, Is.Null,
                "Duplicate spectra warning should be null when no duplicates are present");

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
        public async Task TestKoinaProsit2020IntensityHCDModelFiltersDuplicatePeptidesFromSpectralLibrary()
        {
            var peptides = new List<string> { "PEPTIDEK", "PEPTIDEK", "ELVISLIVESK" };
            var charges = new List<int> { 2, 2, 2 };
            var energies = new List<int> { 35, 35, 35 };
            var retentionTimes = new List<double?> { 100.0, 100.0, 200.0 };
            var model = new Prosit2020IntensityHCD(peptides, charges, energies, retentionTimes, out var warning);
            WarningException? duplicatesWarning = await model.RunInferenceAsync();
            Assert.That(model.PredictedSpectra.Count, Is.EqualTo(2),
                "Duplicate peptide should be filtered out, resulting in 2 unique spectra");
            Assert.That(duplicatesWarning, Is.Not.Null,
                "Duplicate spectra warning should be generated when duplicates are present");
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
            Assert.That(modelHandler.Predictions.Count == numberOfSequences);
        }

        [Test]
        public static async Task TestKoinaProsit2020IntensityHCDModelIsDisposedProperly()
        {
            var peptides = new List<string> { "PEPTIDE", "PEPTIDEK" };
            var charges = new List<int> { 2, 2};
            var energies = new List<int> { 35, 35};
            var retentionTimes = new List<double?> { 100.0, 200.0 };

            // Disposed after use/inference
            var model = new Prosit2020IntensityHCD(peptides, charges, energies, retentionTimes, out var warning);
            Assert.That(warning, Is.Null,
                "Warning should not be generated for valid peptides");
            var duplicatesWarning = await model.RunInferenceAsync();
            Assert.That(duplicatesWarning, Is.Null,
                "Duplicate warning should be null when no duplicates are present");
            Assert.ThrowsAsync<ObjectDisposedException>(async () => await model.RunInferenceAsync(),
                "Running inference on a disposed model should throw ObjectDisposedException");

            // Results should still be accessible after disposal
            Assert.That(model.PredictedSpectra.Count, Is.EqualTo(2),
                "Predicted spectra should still be accessible after model is disposed");

            // Disposed on command
            model = new Prosit2020IntensityHCD(peptides, charges, energies, retentionTimes, out warning);
            Assert.That(warning, Is.Null,
                "Warning should not be generated for valid peptides");
            model.Dispose();
            Assert.ThrowsAsync<ObjectDisposedException>(async () => await model.RunInferenceAsync(),
                "Running inference on a disposed model should throw ObjectDisposedException");
        }

        ///// <summary>
        ///// Reads the Koina input TSV at <c>F:\DiaBenchmark\PXD005573\DiannOut\koina_input.tsv</c>
        ///// (columns: modified_sequence, collision_energy, precursor_charge), obtains iRT predictions
        ///// via Prosit 2019 iRT, then submits the full table for HCD fragment intensity prediction and
        ///// writes a single spectral library to the same directory as the input file.
        /////
        ///// The TSV uses UNIMOD bracket notation (e.g. C[UNIMOD:4], M[UNIMOD:35]); this test converts
        ///// those tokens to mzLib format before submission and converts back for the RT lookup.
        ///// </summary>
        //[Test]
        //public static async Task TestBuildSpectralLibraryFromKoinaInputTsv()
        //{
        //    const string tsvPath = @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.tsv";

        //    if (!File.Exists(tsvPath))
        //        Assert.Ignore("Input TSV not found at the specified path. Update tsvPath to run this test.");

        //    string outDir  = Path.GetDirectoryName(tsvPath)!;
        //    string mspPath = Path.Combine(outDir, Path.GetFileNameWithoutExtension(tsvPath) + ".msp");

        //    // 1. Parse the TSV -------------------------------------------------------
        //    var lines = File.ReadAllLines(tsvPath);
        //    Assert.That(lines.Length, Is.GreaterThan(1), "TSV has no data rows.");

        //    var header    = lines[0].Split('\t');
        //    int seqIdx    = Array.IndexOf(header, "modified_sequence");
        //    int ceIdx     = Array.IndexOf(header, "collision_energy");
        //    int chargeIdx = Array.IndexOf(header, "precursor_charge");

        //    Assert.That(seqIdx    >= 0, "Column 'modified_sequence' not found.");
        //    Assert.That(ceIdx     >= 0, "Column 'collision_energy' not found.");
        //    Assert.That(chargeIdx >= 0, "Column 'precursor_charge' not found.");

        //    // The TSV uses UNIMOD bracket notation written by KoinaTableGenerator.
        //    // mzLib Prosit models require the mzLib modification annotation; convert here.
        //    static string UnimodToMzLib(string seq) => seq
        //        .Replace("[UNIMOD:4]",  "[Common Fixed:Carbamidomethyl on C]")
        //        .Replace("[UNIMOD:35]", "[Common Variable:Oxidation on M]");

        //    var mzLibSequences = new List<string>();
        //    var charges        = new List<int>();
        //    var energies       = new List<int>();

        //    for (int i = 1; i < lines.Length; i++)
        //    {
        //        if (string.IsNullOrWhiteSpace(lines[i])) continue;
        //        var fields = lines[i].Split('\t');
        //        if (fields.Length <= Math.Max(Math.Max(seqIdx, ceIdx), chargeIdx)) continue;

        //        if (!double.TryParse(fields[ceIdx].Trim(),    out double ce))     continue;
        //        if (!int.TryParse(fields[chargeIdx].Trim(),   out int charge))    continue;

        //        mzLibSequences.Add(UnimodToMzLib(fields[seqIdx].Trim()));
        //        charges.Add(charge);
        //        energies.Add((int)ce);
        //    }

        //    Assert.That(mzLibSequences.Count, Is.GreaterThan(0), "No valid rows parsed from TSV.");

        //    // 2. Retention time predictions ------------------------------------------
        //    // Submit unique sequences to avoid redundant API calls.
        //    var uniqueSeqs = mzLibSequences.Distinct().ToList();
        //    var rtModel    = new Prosit2019iRT(uniqueSeqs, out _);
        //    await rtModel.RunInferenceAsync();

        //    // The RT model converts mzLib → UNIMOD internally; predictions are keyed by UNIMOD.
        //    var unimodToRt = new Dictionary<string, double?>(rtModel.Predictions.Count);
        //    for (int i = 0; i < rtModel.Predictions.Count; i++)
        //        unimodToRt[rtModel.Predictions[i].FullSequence] = rtModel.Predictions[i].PredictedRetentionTime;

        //    // Mirror the same conversion to look up RT for each row's mzLib sequence.
        //    static string MzLibToUnimod(string seq) => seq
        //        .Replace("[Common Variable:Oxidation on M]", "[UNIMOD:35]")
        //        .Replace("[Common Fixed:Carbamidomethyl on C]", "[UNIMOD:4]");

        //    var retentionTimes = mzLibSequences
        //        .Select(s => unimodToRt.TryGetValue(MzLibToUnimod(s), out var rt) ? rt : null)
        //        .ToList();

        //    // 3. HCD fragment intensity predictions + write library ------------------
        //    var hcdModel = new Prosit2020IntensityHCD(
        //        mzLibSequences, charges, energies, retentionTimes, out _, mspPath);
        //    await hcdModel.RunInferenceAsync();

        //    // 4. Assertions ----------------------------------------------------------
        //    Assert.That(File.Exists(mspPath),
        //        $"Spectral library was not written to {mspPath}.");
        //    Assert.That(new FileInfo(mspPath).Length, Is.GreaterThan(0),
        //        "Spectral library file is empty.");
        //    Assert.That(hcdModel.PredictedSpectra.Count, Is.GreaterThan(0),
        //        "No spectra were predicted.");
        //    Assert.That(hcdModel.PredictedSpectra.Count, Is.LessThanOrEqualTo(mzLibSequences.Count),
        //        "Cannot have more predicted spectra than input rows.");
        //}
    }
}
