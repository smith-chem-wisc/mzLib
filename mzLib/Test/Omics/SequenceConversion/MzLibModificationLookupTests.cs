using NUnit.Framework;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion;

[TestFixture]
[System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
public class MzLibModificationLookupTests
{
    [Test]
    public void ProteinLookup_DoesNotResolveRnaModifications()
    {
        var lookup = MzLibModificationLookup.ProteinOnly;
        var mod = CanonicalModification.AtResidue(
            residueIndex: 0,
            targetResidue: 'A',
            originalRepresentation: "N6-methyladenosine",
            mzLibId: "N6-methyladenosine");

        var result = lookup.TryResolve(mod);

        Assert.That(result, Is.Null);
    }

    [Test]
    public void RnaLookup_ResolvesRnaModifications()
    {
        var lookup = MzLibModificationLookup.RnaOnly;
        var mod = CanonicalModification.AtResidue(
            residueIndex: 0,
            targetResidue: 'A',
            originalRepresentation: "N6-methyladenosine",
            mzLibId: "N6-methyladenosine");

        var result = lookup.TryResolve(mod);

        Assert.That(result, Is.Not.Null);
        Assert.That(result.Value.IsResolved, Is.True);
        Assert.That(result.Value.MzLibModification, Is.Not.Null);
    }
}
