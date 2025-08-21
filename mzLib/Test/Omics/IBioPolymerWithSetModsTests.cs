using NUnit.Framework;
using Omics;

namespace Test.Omics;

[TestFixture]
public class IBioPolymerWithSetModsTests
{
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
    public void TestGetBaseSequenceFromFullSequence(string fullSequence, string expectedBaseSequence, char startDelimiter, char endDelimiter, char cTerminusDelimiter)
    {
        string actualBaseSequence = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(fullSequence, startDelimiter, endDelimiter, cTerminusDelimiter);
        Assert.That(actualBaseSequence, Is.EqualTo(expectedBaseSequence));
    }
}