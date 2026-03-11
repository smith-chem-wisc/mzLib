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

        var serializer = new UniprotSequenceSerializer();
        var serialized = serializer.Serialize(canonical);

        Assert.That(serialized, Is.EqualTo("PEPS[UniProt:Phosphoserine on S]"));
    }

    [Test]
    public void ThrowsWhenUniProtMappingIsUnavailable()
    {
        var builder = new CanonicalSequenceBuilder("PEP")
            .AddResidueModification(
                residueIndex: 1,
                originalRepresentation: "UnknownMod",
                mass: 42.0);

        var serializer = new UniprotSequenceSerializer();

        Assert.Throws<SequenceConversionException>(() => serializer.Serialize(builder.Build()));
    }
}
