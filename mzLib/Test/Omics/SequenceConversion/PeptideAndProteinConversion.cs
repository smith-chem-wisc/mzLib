using NUnit.Framework;
using Omics;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion;
[TestFixture]
public class PeptideAndProteinConversion
{
    public static IEnumerable<GroundTruthTestData.SequenceConversionTestCase> UniProtTestCases() => GroundTruthTestData.CoreTestCases
        .Concat(GroundTruthTestData.EdgeCases)
        .Where(testCase => !string.IsNullOrWhiteSpace(testCase.UniProtFormat));

    #region PeptideWithSetModifications Conversion Tests

    [Test]
    public static void TestConvertModificationsOnPeptideWithSetModifications()
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
        peptide.ConvertModifications(UniProtModificationLookup.Instance);

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
    public static void TestConvertModificationsOnPeptidePreservesMotifs()
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
        peptide.ConvertModifications(UniProtModificationLookup.Instance);

        // Verify target is preserved
        Assert.That(peptide.AllModsOneIsNterminus[7].Target.ToString(), Is.EqualTo(originalTarget));
        Assert.That(peptide.AllModsOneIsNterminus[7].Target.ToString(), Does.Contain("M"));
    }

    [Test]
    public static void TestConvertModificationsOnPeptideWithCTerminalMod()
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
        peptide.ConvertModifications(UniProtModificationLookup.Instance);

        // Verify C-terminal mod was converted
        Assert.That(peptide.AllModsOneIsNterminus.ContainsKey(9), Is.True);
        Assert.That(peptide.AllModsOneIsNterminus[9].ModificationType, Is.EqualTo("UniProt"));
    }

    [Test]
    public static void TestConvertModificationsOnEmptyPeptide()
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
        Assert.DoesNotThrow(() => peptide.ConvertModifications(UniProtModificationLookup.Instance));

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
        peptide.ConvertModifications(UniProtModificationLookup.Instance);

        var uniprotFormula = peptide.AllModsOneIsNterminus[7].ChemicalFormula;

        // Convert back to MetaMorpheus
        peptide.ConvertModifications(MzLibModificationLookup.Instance);

        var finalFormula = peptide.AllModsOneIsNterminus[7].ChemicalFormula;
        var finalTarget = peptide.AllModsOneIsNterminus[7].Target.ToString();

        // Verify chemistry is preserved
        Assert.That(originalFormula.Equals(uniprotFormula), Is.True);
        Assert.That(originalFormula.Equals(finalFormula), Is.True);
        Assert.That(originalTarget, Is.EqualTo(finalTarget));
    }

    [Test]
    [TestCaseSource(nameof(UniProtTestCases))]
    public static void TestConvertModificationsWithGroundTruthUniProtFormats(GroundTruthTestData.SequenceConversionTestCase testCase)
    {
        var serializer = UniProtSequenceSerializer.Instance;

        var mods = IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(
            testCase.MzLibFormat,
            Mods.AllModsKnownDictionary);

        var protein = new Protein(testCase.ExpectedBaseSequence, "TestProtein");
        var digestionParams = new DigestionParams(protease: "trypsin");
        var peptide = new PeptideWithSetModifications(
            protein,
            digestionParams,
            oneBasedStartResidueInProtein: 1,
            oneBasedEndResidueInProtein: testCase.ExpectedBaseSequence.Length,
            cleavageSpecificity: CleavageSpecificity.Full,
            peptideDescription: "Test",
            missedCleavages: 0,
            allModsOneIsNterminus: mods,
            numFixedMods: 0);

        peptide.ConvertModifications(serializer);

        var canonical = peptide.ToCanonicalSequence();
        var result = serializer.Serialize(canonical, null, serializer.HandlingMode);

        Assert.That(result, Is.EqualTo(testCase.UniProtFormat), testCase.Description);
    }

    #endregion
}
