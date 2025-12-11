using NUnit.Framework;
using Readers.SpectralLibrary;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.IO;
using System.Windows.Media;
using Easy.Common.Extensions;
using System;
using System.ComponentModel;

namespace Test
{
    public class FragmentIntensityPrediction
    {
        [Test]
        public static async Task TestKoinaProsit2020IntensityHCDModelWritesReadableSpectralLibrary()
        {
            var experPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"SpectralLibrary\SpectralLibraryData\myPrositLib.msp");
            var predPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"SpectralLibrary\SpectralLibraryData\koinaTestOutput.msp");
            if (File.Exists(predPath))
            {
                File.Delete(predPath);
            }

            var testLibraryWithoutDecoy = new SpectralLibrary(new List<string> { experPath });
            string aminoacids = @"ACDEFGHIKLMNPQRSTVWY";
            var librarySpectra = testLibraryWithoutDecoy.GetAllLibrarySpectra()
                .Where(p => p.ChargeState < 6 && p.Sequence.Length < 30 && p.Sequence.Length > 1 && p.Sequence.ToHashSet().IsSubsetOf(aminoacids.ToHashSet())).ToList();

            var peptides = librarySpectra.Select(p => p.Sequence).ToList();
            var charges = librarySpectra.Select(p => p.ChargeState).ToList();
            var energies = librarySpectra.Select(p => 35).ToList();
            var retentionTimes = librarySpectra.Select(p => p.RetentionTime).ToList();

            var modelHandler = new Koina.SupportedModels.Prosit2020IntensityHCD.Prosit2020IntensityHCD(peptides, charges, energies, retentionTimes, out var warnings, minIntensityFilter: 1e-6);

            await modelHandler.RunInferenceAsync();
            var predictedSpectra = modelHandler.PredictedSpectra;
            modelHandler.SavePredictedSpectralLibrary(predPath);

            // Test that the predicted spectrum that was saved matches the predicted spectra in memory
            var spectralLibraryTest = new SpectralLibrary(new List<string> { predPath });
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
                    Assert.That(inMemoryFrag.Mz == savedFrag.Mz);
                    // Compare normalized intensities since saved library normalizes to max intensity
                    Assert.That((inMemoryFrag.Intensity / maxInMemoryIntensity) == savedFrag.Intensity);
                    Assert.That(inMemoryFrag.NeutralTheoreticalProduct.ProductType == savedFrag.NeutralTheoreticalProduct.ProductType);
                    Assert.That(inMemoryFrag.NeutralTheoreticalProduct.FragmentNumber == savedFrag.NeutralTheoreticalProduct.FragmentNumber);
                    Assert.That(inMemoryFrag.Charge == savedFrag.Charge);
                }
            }

            testLibraryWithoutDecoy.CloseConnections();
            spectralLibraryTest.CloseConnections();
            File.Delete(predPath);
        }

        [Test]
        public static void TestKoinaProsit2020IntensityHCDModelInputValidation()
        {
            // Arrange
            var validPeptides = new List<string> { "ACDEFGHIK", "LMNPQRSTVWY" };
            var validCharges = new List<int> { 2, 2 };
            var validEnergies = new List<int> { 35, 35 };
            var validRetentionTimes = new List<double?> { 100.0, 200.0 };
            WarningException warning;

            // Mismatched lengths
            var mismatchedCharges = new List<int> { 2 };
            Assert.Throws<System.ArgumentException>(() =>
                new Koina.SupportedModels.Prosit2020IntensityHCD.Prosit2020IntensityHCD(validPeptides, mismatchedCharges, validEnergies, validRetentionTimes, out var _));

            // Invalid charge state
            var invalidCharges = new List<int> { 7, 2 };
            var model = new Koina.SupportedModels.Prosit2020IntensityHCD.Prosit2020IntensityHCD(validPeptides, invalidCharges, validEnergies, validRetentionTimes, out warning);
            Assert.That(warning, Is.Not.Null);
            Assert.That(model.PeptideSequences.Count == 1);

            // Invalid peptide sequences (too long, too short)
            var invalidPeptides = new List<string> { "PEPTIDAAAAAAAAAAAAAAAAAAAAAAAAA", "" };
            model = new Koina.SupportedModels.Prosit2020IntensityHCD.Prosit2020IntensityHCD(invalidPeptides, validCharges, validEnergies, validRetentionTimes, out warning);
            Assert.That(warning, Is.Not.Null);
            Assert.That(model.PeptideSequences.Count == 0);

            // Valid peptide with modifications (length > 30 with mods, but valid without mods)
            var peptidesWithMods = new List<string> { "M[Common Variable:Oxidation on M]PEPTIDEC[Common Fixed:Carbamidomethyl on C]A", "LMNPQRSTVWY" };
            model = new Koina.SupportedModels.Prosit2020IntensityHCD.Prosit2020IntensityHCD(peptidesWithMods, validCharges, validEnergies, validRetentionTimes, out warning);
            Assert.That(warning, Is.Null);
            Assert.That(model.PeptideSequences.Count == 2);
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
            var modelHandler = new Koina.SupportedModels.Prosit2020IntensityHCD.Prosit2020IntensityHCD(peptides, charges, energies, retentionTimes, out var warnings);
            Assert.DoesNotThrowAsync(async () => await modelHandler.RunInferenceAsync());
            Assert.That(modelHandler.PredictedSpectra.Count == numberOfSequences);
        }
    }
}
