using NUnit.Framework;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion;

[TestFixture]
public class CanonicalSequenceTests
{
    [Test]
    public void WithModification_SortsByPositionTypeAndResidueIndex()
    {
        var sequence = CanonicalSequence.Unmodified("PEPTIDE")
            .WithModification(CanonicalModification.AtResidue(3, 'T', "r3", mass: 1.0))
            .WithModification(CanonicalModification.AtCTerminus("c", mass: 2.0))
            .WithModification(CanonicalModification.AtResidue(1, 'E', "r1", mass: 3.0))
            .WithModification(CanonicalModification.AtNTerminus("n", mass: 4.0));

        Assert.That(sequence.Modifications[0].PositionType, Is.EqualTo(ModificationPositionType.NTerminus));
        Assert.That(sequence.Modifications[1].ResidueIndex, Is.EqualTo(1));
        Assert.That(sequence.Modifications[2].ResidueIndex, Is.EqualTo(3));
        Assert.That(sequence.Modifications[3].PositionType, Is.EqualTo(ModificationPositionType.CTerminus));
    }

    [Test]
    public void TotalModificationMass_ReturnsNullWhenAnyMassMissing()
    {
        var withMissingMass = CanonicalSequence.Unmodified("PEPTIDE")
            .WithModification(CanonicalModification.AtResidue(2, 'P', "known", mass: 10.0))
            .WithModification(CanonicalModification.AtResidue(4, 'I', "unknown"));

        Assert.That(withMissingMass.TotalModificationMass, Is.Null);
        Assert.That(withMissingMass.AllModificationsHaveMass, Is.False);
    }

    [Test]
    public void ToMarkedSequence_HandlesTerminalAndResidueMods()
    {
        var sequence = CanonicalSequence.Unmodified("PEP")
            .WithModification(CanonicalModification.AtNTerminus("n", mass: 1.0))
            .WithModification(CanonicalModification.AtResidue(1, 'E', "r", mass: 2.0))
            .WithModification(CanonicalModification.AtCTerminus("c", mass: 3.0));

        Assert.That(sequence.ToMarkedSequence(), Is.EqualTo("[*]-PE[*]P-[*]"));
    }

    [Test]
    public void GetModificationAt_ReturnsExpectedModification()
    {
        var sequence = CanonicalSequence.Unmodified("PEPTIDE")
            .WithModification(CanonicalModification.AtResidue(2, 'P', "mod", mass: 15.99));

        Assert.That(sequence.HasModificationAt(2), Is.True);
        Assert.That(sequence.GetModificationAt(2), Is.Not.Null);
        Assert.That(sequence.GetModificationAt(0), Is.Null);
        Assert.That(sequence.TotalModificationMass, Is.EqualTo(15.99).Within(1e-9));
    }
}
