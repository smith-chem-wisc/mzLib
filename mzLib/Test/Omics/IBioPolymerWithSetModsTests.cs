using NUnit.Framework;
using Omics;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using System.Diagnostics.CodeAnalysis;

namespace Test.Omics;

/// <summary>
/// Tests for IBioPolymerWithSetMods static utility methods.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class IBioPolymerWithSetModsTests
{
    /// <summary>
    /// Verifies GetBaseSequenceFromFullSequence correctly strips modification annotations from sequences.
    /// Critical: Base sequence extraction is required for peptide matching, protein inference,
    /// and sequence coverage calculations. Supports various delimiter formats used by different tools.
    /// </summary>
    /// <param name="fullSequence">Modified sequence with annotations</param>
    /// <param name="expectedBaseSequence">Expected sequence without modifications</param>
    /// <param name="startDelimiter">Opening delimiter for modifications (e.g., '[', '{', '&lt;', '(')</param>
    /// <param name="endDelimiter">Closing delimiter for modifications (e.g., ']', '}', '&gt;', ')')</param>
    /// <param name="cTerminusDelimiter">Delimiter separating C-terminal modifications (e.g., '-')</param>
    [Test]
    [TestCase("AC[Phospho]DEFGH", "ACDEFGH", '[', ']', '-')]
    [TestCase("A[Oxidation]C[Phospho]DEFGH", "ACDEFGH", '[', ']', '-')]
    [TestCase("A[Oxidation]C[Phospho]D[Carbamidomethyl]EFGH", "ACDEFGH", '[', ']', '-')]
    [TestCase("A[Oxidation]C[Phospho]D[Carbamidomethyl]E[Acetyl]FGH", "ACDEFGH", '[', ']', '-')]
    [TestCase("A[Oxidation]C[Phospho]D[Carbamidomethyl]E[Acetyl]F[Phospho]GH", "ACDEFGH", '[', ']', '-')]
    [TestCase("A[Oxidation]C[Phospho]D[Carbamidomethyl]E[Acetyl]F[Phospho]G[Oxidation]H", "ACDEFGH", '[', ']', '-')]
    [TestCase("AC{Phospho}DEFGH", "ACDEFGH", '{', '}', '-')]
    [TestCase("A<Oxidation>C<Phospho>DEFGH", "ACDEFGH", '<', '>', '-')]
    [TestCase("A(Oxidation)C(Phospho)DEFGH", "ACDEFGH", '(', ')', '-')]
    [TestCase("ACDEFGH-[TestMod: ModNameX]", "ACDEFGH", '[', ']', '-')]
    [TestCase("ACDEFGH&[TestMod: ModNameX]", "ACDEFGH", '[', ']', '&')]
    [TestCase("[TestMod: ModName]ACD[TestMod1: ModNameX]EFGH[TestMod2: ModNameY]-[TestMod3: ModNameZ]", "ACDEFGH", '[', ']', '-')]
    public void GetBaseSequenceFromFullSequence_StripsModifications(
        string fullSequence,
        string expectedBaseSequence,
        char startDelimiter,
        char endDelimiter,
        char cTerminusDelimiter)
    {
        string actualBaseSequence = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(
            fullSequence, startDelimiter, endDelimiter, cTerminusDelimiter);

        Assert.That(actualBaseSequence, Is.EqualTo(expectedBaseSequence));
    }

    [Test]
    public void GetModificationDictionaryFromFullSequence_ResolvesUnimodModification()
    {
        var mods = IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(
            "PEPC[UNIMOD:4]IDE",
            Mods.AllKnownProteinModsDictionary);

        Assert.That(mods, Contains.Key(5));
        Assert.That(mods[5].IdWithMotif, Is.EqualTo("Carbamidomethyl on C"));
    }

    [Test]
    public void GetModificationDictionaryFromFullSequence_UsesResidueToDisambiguateUnimod()
    {
        var mods = IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(
            "PEPT[UNIMOD:21]IDE",
            Mods.AllKnownProteinModsDictionary);

        Assert.That(mods, Contains.Key(5));
        Assert.That(mods[5].Target?.ToString(), Does.Contain("T"));
    }

    [Test]
    public void PeptideWithSetModifications_Ctor_AcceptsUnimodSequence()
    {
        var peptide = new PeptideWithSetModifications("PEPC[UNIMOD:4]IDE");

        Assert.That(peptide.BaseSequence, Is.EqualTo("PEPCIDE"));
        Assert.That(peptide.AllModsOneIsNterminus, Contains.Key(5));
        Assert.That(peptide.AllModsOneIsNterminus[5].IdWithMotif, Is.EqualTo("Carbamidomethyl on C"));
    }

    [Test]
    [TestCase("-XYZ--ABC")]
    [TestCase("AB--------")]
    [TestCase("-----F----*")]
    public void GetModificationDictionaryFromFullSequence_LegacyDashPlaceholders_ReturnsEmptyDictionary(string fullSequence)
    {
        var mods = IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(
            fullSequence,
            Mods.AllKnownProteinModsDictionary);

        Assert.That(mods, Is.Empty);
    }

    [Test]
    [TestCase("EKVLTSSAR(2)")]
    [TestCase("EKVLTSSAR(2)SLGKVGTR(4)")]
    [TestCase("LLDNAAADLAAISGQKPLITKAR(21)")]
    public void GetModificationDictionaryFromFullSequence_CrosslinkAnnotatedSequence_ReturnsEmptyDictionary(string fullSequence)
    {
        var mods = IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(
            fullSequence,
            Mods.AllKnownProteinModsDictionary);

        Assert.That(mods, Is.Empty);
    }
}
