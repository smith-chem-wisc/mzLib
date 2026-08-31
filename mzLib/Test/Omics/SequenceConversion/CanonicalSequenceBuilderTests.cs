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

    [Test]
    public void AddNTerminalModification_UsesFirstResidueForTarget()
    {
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddNTerminalModification("Acetyl")
            .Build();

        Assert.That(canonical.NTerminalModification.HasValue, Is.True);
        Assert.That(canonical.NTerminalModification!.Value.TargetResidue, Is.EqualTo('P'));
    }

    [Test]
    public void AddCTerminalModification_UsesLastResidueForTarget()
    {
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddCTerminalModification("Amidate")
            .Build();

        Assert.That(canonical.CTerminalModification.HasValue, Is.True);
        Assert.That(canonical.CTerminalModification!.Value.TargetResidue, Is.EqualTo('E'));
    }

    [Test]
    public void Reset_ClearsBuilderState()
    {
        var builder = new CanonicalSequenceBuilder("PEP")
            .WithSourceFormat("Custom")
            .AddNTerminalModification("Acetyl");

        builder.Reset();

        Assert.That(builder.BaseSequence, Is.EqualTo(string.Empty));
        Assert.That(builder.ModificationCount, Is.Zero);
        Assert.That(builder.SourceFormat, Is.EqualTo("Unknown"));
    }

    [Test]
    public void AddResidueModification_InsertsResidueModification()
    {
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(2, "Oxidation")
            .Build();

        var mod = canonical.GetModificationAt(2);

        Assert.That(mod.HasValue, Is.True);
        Assert.That(mod!.Value.OriginalRepresentation, Is.EqualTo("Oxidation"));
        Assert.That(mod.Value.TargetResidue, Is.EqualTo('P'));
        Assert.That(canonical.ModificationCount, Is.EqualTo(1));
    }

    [Test]
    public void AddModifications_AppendsEachModification()
    {
        var additionalMods = new[]
        {
            CanonicalModification.AtResidue(1, 'E', "Label"),
            CanonicalModification.AtCTerminus("Cterm")
        };

        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddModifications(additionalMods)
            .Build();

        Assert.That(canonical.ModificationCount, Is.EqualTo(2));
        Assert.That(canonical.GetModificationAt(1).Value.OriginalRepresentation, Is.EqualTo("Label"));
        Assert.That(canonical.CTerminalModification.HasValue, Is.True);
        Assert.That(canonical.CTerminalModification!.Value.OriginalRepresentation, Is.EqualTo("Cterm"));
    }

    [Test]
    public void RemoveModificationAt_RemovesOnlyMatchingIndex()
    {
        var builder = new CanonicalSequenceBuilder("PEP")
            .AddResidueModification(0, "First")
            .AddResidueModification(1, "Second");

        builder.RemoveModificationAt(0);

        var canonical = builder.Build();

        Assert.That(canonical.ModificationCount, Is.EqualTo(1));
        Assert.That(canonical.GetModificationAt(0).HasValue, Is.False);
        Assert.That(canonical.GetModificationAt(1).Value.OriginalRepresentation, Is.EqualTo("Second"));

        builder.RemoveModificationAt(5); // out of range
        Assert.That(builder.ModificationCount, Is.EqualTo(1));
    }

    [Test]
    public void RemoveModificationsAtResidue_ClearsResidueModifications()
    {
        var builder = new CanonicalSequenceBuilder("PEP")
            .AddResidueModification(0, "First")
            .AddResidueModification(0, "Alternate")
            .AddResidueModification(2, "Last");

        builder.RemoveModificationsAtResidue(0);

        var canonical = builder.Build();

        Assert.That(canonical.GetModificationAt(0).HasValue, Is.False);
        Assert.That(canonical.ModificationCount, Is.EqualTo(1));
        Assert.That(canonical.GetModificationAt(2).Value.OriginalRepresentation, Is.EqualTo("Last"));
    }

    [Test]
    public void ClearModifications_PreservesSequence()
    {
        var builder = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(0, "First")
            .AddCTerminalModification("Cterm");

        var canonical = builder.ClearModifications().Build();

        Assert.That(canonical.ModificationCount, Is.Zero);
        Assert.That(canonical.BaseSequence, Is.EqualTo("PEPTIDE"));
    }

    [Test]
    public void Constructor_CopiesExistingSequenceAndSourceFormat()
    {
        var original = new CanonicalSequenceBuilder("PEP")
            .AddNTerminalModification("Acetyl")
            .AddResidueModification(1, "Label")
            .WithSourceFormat("MzLib")
            .Build();

        var builder = new CanonicalSequenceBuilder(original);
        var duplicate = builder.AddResidueModification(2, "Extra").Build();

        Assert.That(duplicate.BaseSequence, Is.EqualTo(original.BaseSequence));
        Assert.That(duplicate.SourceFormat, Is.EqualTo("MzLib"));
        Assert.That(duplicate.ModificationCount, Is.EqualTo(original.ModificationCount + 1));
        Assert.That(original.ModificationCount, Is.EqualTo(2));
    }

    [Test]
    public void SourceFormatSetter_NormalizesNull()
    {
        var builder = new CanonicalSequenceBuilder("PEP")
        {
            SourceFormat = "Custom"
        };

        builder.SourceFormat = null!;

        Assert.That(builder.SourceFormat, Is.EqualTo("Unknown"));
        Assert.That(builder.Build().SourceFormat, Is.EqualTo("Unknown"));
    }

    [Test]
    public void Validate_ReturnsErrorForEmptyBaseSequence()
    {
        var builder = new CanonicalSequenceBuilder();

        var errors = builder.Validate();

        Assert.That(errors, Has.Some.EqualTo("Base sequence is empty."));
        Assert.That(builder.IsValid, Is.False);
    }
}
