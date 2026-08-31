using NUnit.Framework;
using Omics;
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
}