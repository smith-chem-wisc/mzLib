using System.Linq;
using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion;

[TestFixture]
[System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
public class UniprotSequenceSerializerTests
{
    [Test]
    public void SerializesMetaMorpheusModificationsToUniProtNotation()
    {
        var metaMod = Mods.AllKnownProteinModsDictionary["Phosphorylation on S"];

        var builder = new CanonicalSequenceBuilder("PEPS")
            .AddResidueModification(
                residueIndex: 3,
                originalRepresentation: "Common Biological:Phosphorylation on S",
                mass: metaMod.MonoisotopicMass,
                formula: metaMod.ChemicalFormula,
                mzLibId: metaMod.IdWithMotif,
                mzLibModification: metaMod);

        var canonical = builder.Build();

        var serializer = new UniProtSequenceSerializer();
        var serialized = serializer.Serialize(canonical);

        Assert.That(serialized, Is.EqualTo("PEPS[Phosphoserine]"));
    }

    [Test]
    public void ThrowsWhenUniProtMappingIsUnavailable()
    {
        var builder = new CanonicalSequenceBuilder("PEP")
            .AddResidueModification(
                residueIndex: 1,
                originalRepresentation: "UnknownMod",
                mass: 42.0);

        var serializer = new UniProtSequenceSerializer();

        Assert.Throws<SequenceConversionException>(() => serializer.Serialize(builder.Build()));
    }

    [Test]
    public void SerializesMetaMorpheusThreonineToUniProtNotation()
    {
        var metaMod = Mods.AllKnownProteinModsDictionary["Phosphorylation on T"];

        var builder = new CanonicalSequenceBuilder("PEPT")
            .AddResidueModification(
                residueIndex: 3,
                originalRepresentation: "Common Biological:Phosphorylation on T",
                mass: metaMod.MonoisotopicMass,
                formula: metaMod.ChemicalFormula,
                mzLibId: metaMod.IdWithMotif,
                mzLibModification: metaMod);

        var serializer = new UniProtSequenceSerializer();
        var serialized = serializer.Serialize(builder.Build());

        Assert.That(serialized, Is.EqualTo("PEPT[Phosphothreonine]"));
    }

    [Test]
    public void Lookup_TrimsSuffixesAndResolvesUniProtName()
    {
        var lookup = new UniProtModificationLookup();
        var canonical = CanonicalModification.AtResidue(
            residueIndex: 0,
            targetResidue: 'S',
            originalRepresentation: "UniProt:Phosphoserine-L- on S");

        var resolved = lookup.TryResolve(canonical);

        Assert.That(resolved, Is.Not.Null);
        Assert.That(resolved!.Value.MzLibModification, Is.Not.Null);
        Assert.That(resolved.Value.MzLibModification!.IdWithMotif, Does.Contain("Phosphoserine"));
    }

    [Test]
    public void Serializer_RemovesUnknownModificationsWhenConfigured()
    {
        var builder = new CanonicalSequenceBuilder("PEP")
            .AddResidueModification(
                residueIndex: 1,
                originalRepresentation: "UnknownMod",
                mass: 15.0);

        var serializer = new UniProtSequenceSerializer();
        var warnings = new ConversionWarnings();

        var result = serializer.Serialize(builder.Build(), warnings, SequenceConversionHandlingMode.RemoveIncompatibleElements);

        Assert.That(result, Is.EqualTo("PEP"));
        Assert.That(warnings.IncompatibleItems, Is.Not.Empty);
    }
}
