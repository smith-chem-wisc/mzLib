using NUnit.Framework;
using Omics.SequenceConversion;
using System;
using System.Linq;

namespace Test.Omics.SequenceConversion;

[TestFixture]
public class CanonicalSequenceBuilderTests
{
    [Test]
    public void AddResidueModification_InvalidIndex_Throws()
    {
        var builder = new CanonicalSequenceBuilder("PEP");

        Assert.That(() => builder.AddResidueModification(5, "Mod"), Throws.InstanceOf<ArgumentOutOfRangeException>());
    }

    [Test]
    public void Build_SortsModificationsByPosition()
    {
        var builder = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(3, "Third")
            .AddCTerminalModification("CTerm")
            .AddResidueModification(0, "First")
            .AddNTerminalModification("NTerm");

        var canonical = builder.Build();
        var mods = canonical.Modifications.ToArray();

        Assert.That(mods[0].PositionType, Is.EqualTo(ModificationPositionType.NTerminus));
        Assert.That(mods[1].ResidueIndex, Is.EqualTo(0));
        Assert.That(mods[2].ResidueIndex, Is.EqualTo(3));
        Assert.That(mods[^1].PositionType, Is.EqualTo(ModificationPositionType.CTerminus));
    }

    [Test]
    public void Validate_DetectsMultipleErrorConditions()
    {
        var builder = new CanonicalSequenceBuilder("PEP");
        builder.AddModification(CanonicalModification.AtResidue(5, 'K', "OutOfRange"));
        builder.AddModification(CanonicalModification.AtResidue(1, 'K', "TargetMismatch"));
        builder.AddModification(CanonicalModification.AtResidue(1, 'E', "DuplicateOne"));
        builder.AddModification(CanonicalModification.AtResidue(1, 'E', "DuplicateTwo"));
        builder.AddModification(CanonicalModification.AtNTerminus("NTermOne"));
        builder.AddModification(CanonicalModification.AtNTerminus("NTermTwo"));
        builder.AddModification(CanonicalModification.AtCTerminus("CTermOne"));
        builder.AddModification(CanonicalModification.AtCTerminus("CTermTwo"));

        var errors = builder.Validate();

        Assert.That(errors, Has.Some.Contains("invalid index"));
        Assert.That(errors, Has.Some.Contains("targets 'K'"));
        Assert.That(errors, Has.Some.Contains("Multiple modifications at residue index 1"));
        Assert.That(errors, Has.Some.Contains("Multiple N-terminal modifications"));
        Assert.That(errors, Has.Some.Contains("Multiple C-terminal modifications"));
    }
}
