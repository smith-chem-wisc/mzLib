using NUnit.Framework;
using Readers;
using Omics.Modifications;
using Chemistry;
using System.Linq;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using Omics.Digestion;

namespace Test.FileReadingTests
{
    public static class ModificationConversions
    {
        [Test]
        public static void ModsLoadWithoutCrash()
        {
            var allMods = ModificationConverter.AllKnownMods;
            Assert.That(allMods.Count, Is.GreaterThanOrEqualTo(3110)); // number at time of creation
        }

        [Test]
        [TestCase("Carbamidomethyl", "Carbamidomethyl", 'C')]
        [TestCase("O-(ADP-ribosyl)-L-serine", "ADP-ribosylation", 'S')]
        [TestCase("Acetyl", "Acetylation", 'K')]
        [TestCase("N6,N6,N6-trimethyl-L-lysine", "N6,N6,N6-trimethyllysine", 'K')]
        [TestCase("N,N,N-trimethyl-L-alanine", "N,N,N-trimethylalanine", 'A')]
        [TestCase("O4'-(phospho-5'-adenosine)-L-tyrosine", "Phosphoadenosine", 'Y')]
        [TestCase("O4'-(phopho-ChickenPotPie-adenosine)-L-tyrosine", "Phosphoadenosine", 'Y')]
        [TestCase("Acetyldeoxyhypusine", "Acetyldeoxyhypusine", 'K')]
        [TestCase("15N-oxobutanoic", "15N-oxobutanoic", 'T')]
        [TestCase("Amidine", "Amidine", 'T')]
        public static void ModsConvertToModification(string name, string expectedId, char residue)
        {
            var modification = ModificationConverter.GetClosestMod(name, residue, null);
            Assert.That(modification.OriginalId, Is.EqualTo(expectedId));

            if (name == "Amidine") // Edge Case where we want to pass in the wrong motif and check if it got the right one
                residue = 'X';

            var withMotif = $"{expectedId} on {residue}";
            Assert.That(modification.IdWithMotif, Does.Contain(withMotif));
        }

        #region PeptideWithSetModifications Conversion Tests

        [Test]
        public static void TestConvertModsOnPeptideWithSetModifications()
        {
            // Create a protein with MetaMorpheus-style modifications
            ModificationMotif.TryGetMotif("S", out var motifS);
            ModificationMotif.TryGetMotif("K", out var motifK);
            
            var phosphoMM = new Modification(
                _originalId: "Phosphorylation",
                _modificationType: "Common Biological",
                _target: motifS,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));
            
            var acetylMM = new Modification(
                _originalId: "Acetylation",
                _modificationType: "Common Biological", 
                _target: motifK,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"));

            var protein = new Protein("PEPTKSDE", "TestProtein");
            
            var modsOneIsNterm = new Dictionary<int, Modification>
            {
                { 1, acetylMM }, // N-terminal acetylation
                { 4, acetylMM }, // Acetyl on K at position 4 (P-E-P-T-K)
                { 6, phosphoMM }  // Phospho on S at position 6 (P-E-P-T-K-S)
            };

            var digestionParams = new DigestionParams(protease: "trypsin");
            var peptide = new PeptideWithSetModifications(
                protein, 
                digestionParams, 
                oneBasedStartResidueInProtein: 1,
                oneBasedEndResidueInProtein: 8,
                cleavageSpecificity: CleavageSpecificity.Full,
                peptideDescription: "Test",
                missedCleavages: 0,
                allModsOneIsNterminus: modsOneIsNterm,
                numFixedMods: 0);

            // Verify original mods are MetaMorpheus style
            Assert.That(peptide.AllModsOneIsNterminus[1].ModificationType, Is.EqualTo("Common Biological"));
            Assert.That(peptide.AllModsOneIsNterminus[4].ModificationType, Is.EqualTo("Common Biological"));
            Assert.That(peptide.AllModsOneIsNterminus[6].ModificationType, Is.EqualTo("Common Biological"));

            // Convert to UniProt convention
            peptide.ConvertMods(ModificationNamingConvention.UniProt);

            // Verify conversions
            Assert.That(peptide.AllModsOneIsNterminus[1].ModificationType, Is.EqualTo("UniProt"));
            Assert.That(peptide.AllModsOneIsNterminus[4].ModificationType, Is.EqualTo("UniProt"));
            Assert.That(peptide.AllModsOneIsNterminus[6].ModificationType, Is.EqualTo("UniProt"));

            // Verify chemical formulas are preserved
            Assert.That(peptide.AllModsOneIsNterminus[1].ChemicalFormula, Is.Not.Null);
            Assert.That(peptide.AllModsOneIsNterminus[4].ChemicalFormula, Is.Not.Null);
            Assert.That(peptide.AllModsOneIsNterminus[6].ChemicalFormula, Is.Not.Null);
        }

        [Test]
        public static void TestConvertModsOnPeptidePreservesMotifs()
        {
            // Test that conversion preserves amino acid targets
            ModificationMotif.TryGetMotif("M", out var motifM);
            
            var oxidationMM = new Modification(
                _originalId: "Oxidation",
                _modificationType: "Common Variable",
                _target: motifM,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: ChemicalFormula.ParseFormula("O1"));

            var protein = new Protein("PEPTMIDE", "TestProtein");
            
            var modsOneIsNterm = new Dictionary<int, Modification>
            {
                { 7, oxidationMM } // Oxidation on M
            };

            var digestionParams = new DigestionParams(protease: "trypsin");
            var peptide = new PeptideWithSetModifications(
                protein, 
                digestionParams, 
                oneBasedStartResidueInProtein: 1,
                oneBasedEndResidueInProtein: 8,
                cleavageSpecificity: CleavageSpecificity.Full,
                peptideDescription: "Test",
                missedCleavages: 0,
                allModsOneIsNterminus: modsOneIsNterm,
                numFixedMods: 0);

            // Get original target
            var originalTarget = peptide.AllModsOneIsNterminus[7].Target.ToString();
            
            // Convert to UniProt
            peptide.ConvertMods(ModificationNamingConvention.UniProt);

            // Verify target is preserved
            Assert.That(peptide.AllModsOneIsNterminus[7].Target.ToString(), Is.EqualTo(originalTarget));
            Assert.That(peptide.AllModsOneIsNterminus[7].Target.ToString(), Does.Contain("M"));
        }

        [Test]
        public static void TestConvertModsOnPeptideWithCTerminalMod()
        {
            // Test conversion with C-terminal modification
            ModificationMotif.TryGetMotif("X", out var motifX);
            
            var amidationMM = new Modification(
                _originalId: "Amidation",
                _modificationType: "Common Biological",
                _target: motifX,
                _locationRestriction: "Peptide C-terminal.",
                _chemicalFormula: ChemicalFormula.ParseFormula("H1N1"));

            var protein = new Protein("PEPTIDE", "TestProtein");
            
            var modsOneIsNterm = new Dictionary<int, Modification>
            {
                { 9, amidationMM } // C-terminal mod (length 7 + 2)
            };

            var digestionParams = new DigestionParams(protease: "trypsin");
            var peptide = new PeptideWithSetModifications(
                protein, 
                digestionParams, 
                oneBasedStartResidueInProtein: 1,
                oneBasedEndResidueInProtein: 7,
                cleavageSpecificity: CleavageSpecificity.Full,
                peptideDescription: "Test",
                missedCleavages: 0,
                allModsOneIsNterminus: modsOneIsNterm,
                numFixedMods: 0);

            // Convert to UniProt
            peptide.ConvertMods(ModificationNamingConvention.UniProt);

            // Verify C-terminal mod was converted
            Assert.That(peptide.AllModsOneIsNterminus.ContainsKey(9), Is.True);
            Assert.That(peptide.AllModsOneIsNterminus[9].ModificationType, Is.EqualTo("UniProt"));
        }

        [Test]
        public static void TestConvertModsOnEmptyPeptide()
        {
            // Test that conversion works with no modifications
            var protein = new Protein("PEPTIDE", "TestProtein");
            
            var modsOneIsNterm = new Dictionary<int, Modification>();

            var digestionParams = new DigestionParams(protease: "trypsin");
            var peptide = new PeptideWithSetModifications(
                protein, 
                digestionParams, 
                oneBasedStartResidueInProtein: 1,
                oneBasedEndResidueInProtein: 7,
                cleavageSpecificity: CleavageSpecificity.Full,
                peptideDescription: "Test",
                missedCleavages: 0,
                allModsOneIsNterminus: modsOneIsNterm,
                numFixedMods: 0);

            // Should not throw
            Assert.DoesNotThrow(() => peptide.ConvertMods(ModificationNamingConvention.UniProt));
            
            // Should still have no mods
            Assert.That(peptide.AllModsOneIsNterminus.Count, Is.EqualTo(0));
        }

        [Test]
        public static void TestPeptideConversionRoundTrip()
        {
            // Test converting from MetaMorpheus to UniProt and back preserves chemistry
            ModificationMotif.TryGetMotif("S", out var motifS);
            
            var phosphoMM = new Modification(
                _originalId: "Phosphorylation",
                _modificationType: "Common Biological",
                _target: motifS,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));

            var protein = new Protein("PEPTSIDE", "TestProtein");
            
            var modsOneIsNterm = new Dictionary<int, Modification>
            {
                { 7, phosphoMM }
            };

            var digestionParams = new DigestionParams(protease: "trypsin");
            var peptide = new PeptideWithSetModifications(
                protein, 
                digestionParams, 
                oneBasedStartResidueInProtein: 1,
                oneBasedEndResidueInProtein: 8,
                cleavageSpecificity: CleavageSpecificity.Full,
                peptideDescription: "Test",
                missedCleavages: 0,
                allModsOneIsNterminus: modsOneIsNterm,
                numFixedMods: 0);

            var originalFormula = peptide.AllModsOneIsNterminus[7].ChemicalFormula;
            var originalTarget = peptide.AllModsOneIsNterminus[7].Target.ToString();

            // Convert to UniProt
            peptide.ConvertMods(ModificationNamingConvention.UniProt);
            
            var uniprotFormula = peptide.AllModsOneIsNterminus[7].ChemicalFormula;
            
            // Convert back to MetaMorpheus
            peptide.ConvertMods(ModificationNamingConvention.MetaMorpheus);
            
            var finalFormula = peptide.AllModsOneIsNterminus[7].ChemicalFormula;
            var finalTarget = peptide.AllModsOneIsNterminus[7].Target.ToString();

            // Verify chemistry is preserved
            Assert.That(originalFormula.Equals(uniprotFormula), Is.True);
            Assert.That(originalFormula.Equals(finalFormula), Is.True);
            Assert.That(originalTarget, Is.EqualTo(finalTarget));
        }

        #endregion
    }
}
