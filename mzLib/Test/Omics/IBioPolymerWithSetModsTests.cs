using NUnit.Framework;
using Omics;

namespace Test.Omics;

[TestFixture]
public class IBioPolymerWithSetModsTests
{
    [Test]
    [TestCase("AC[Phospho]DEFGH", "ACDEFGH", '[', ']')]
    [TestCase("A[Oxidation]C[Phospho]DEFGH", "ACDEFGH", '[', ']')]
    [TestCase("A[Oxidation]C[Phospho]D[Carbamidomethyl]EFGH", "ACDEFGH", '[', ']')]
    [TestCase("A[Oxidation]C[Phospho]D[Carbamidomethyl]E[Acetyl]FGH", "ACDEFGH", '[', ']')]
    [TestCase("A[Oxidation]C[Phospho]D[Carbamidomethyl]E[Acetyl]F[Phospho]GH", "ACDEFGH", '[', ']')]
    [TestCase("A[Oxidation]C[Phospho]D[Carbamidomethyl]E[Acetyl]F[Phospho]G[Oxidation]H", "ACDEFGH", '[', ']')]
    [TestCase("AC{Phospho}DEFGH", "ACDEFGH", '{', '}')]
    [TestCase("A<Oxidation>C<Phospho>DEFGH", "ACDEFGH", '<', '>')]
    [TestCase("A(Oxidation)C(Phospho)DEFGH", "ACDEFGH", '(', ')')]
    public void TestGetBaseSequenceFromFullSequence(string fullSequence, string expectedBaseSequence, char startDelimiter, char endDelimiter)
    {
        string actualBaseSequence = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(fullSequence, startDelimiter, endDelimiter);
        Assert.That(actualBaseSequence, Is.EqualTo(expectedBaseSequence));
    }
}