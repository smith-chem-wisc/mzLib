using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Omics.Modifications;
using Omics.Modifications.Conversion;

namespace Test.Omics;

[TestFixture]
[ExcludeFromCodeCoverage]
public class ModificationCrossRefIndexTests
{
    [Test]
    public void ModificationCrossRefIndex_GetByDatabaseId_ReturnsMatchingMod()
    {
        var modWithResid = Mods.AllKnownMods.FirstOrDefault(mod =>
            mod.DatabaseReference != null
            && mod.DatabaseReference.TryGetValue("RESID", out var ids)
            && ids.Any());

        Assert.That(modWithResid, Is.Not.Null, "Expected at least one modification with RESID reference.");

        var resid = modWithResid!.DatabaseReference!["RESID"].First();
        var matches = ModificationCrossRefIndex.Global.GetByDatabaseId("RESID", resid);

        Assert.That(matches, Does.Contain(modWithResid));
    }

    [Test]
    public void ModificationCrossRefIndex_GetByAccession_ReturnsMatchingMod()
    {
        var modWithAccession = Mods.AllKnownMods.FirstOrDefault(mod => !string.IsNullOrWhiteSpace(mod.Accession));
        Assert.That(modWithAccession, Is.Not.Null, "Expected at least one modification with an accession.");

        var matches = ModificationCrossRefIndex.Global.GetByAccession(modWithAccession!.Accession);

        Assert.That(matches, Does.Contain(modWithAccession));
    }
}
