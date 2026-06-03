using System.Linq;
using NUnit.Framework;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion;

[TestFixture]
public class ChronologerSequenceSerializerTests
{
    private ChronologerSequenceSerializer _serializer = null!;

    [SetUp]
    public void SetUp()
    {
        _serializer = new ChronologerSequenceSerializer();
    }

    [Test]
    public void SchemaReferencesSingleton()
    {
        Assert.That(_serializer.Schema, Is.SameAs(ChronologerSequenceFormatSchema.Instance));
        Assert.That(_serializer.Schema.FormatName, Is.EqualTo(_serializer.FormatName));
    }

    [Test]
    public void CanSerialize_RejectsNonCanonicalAminoAcids()
    {
        var canonical = CanonicalSequence.Unmodified("AX");

        Assert.That(_serializer.CanSerialize(canonical), Is.False);
    }

    [Test]
    public void Serialize_ReturnNullForEmptySequenceAndRecordsFailure()
    {
        var warnings = new ConversionWarnings();

        var result = _serializer.Serialize(CanonicalSequence.Empty, warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(result, Is.Null);
        Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.InvalidSequence));
    }

    [Test]
    public void Serialize_ThrowsWhenSequenceExceedsMaxLength()
    {
        var longSequence = new string('A', ChronologerSequenceFormatSchema.MaxSequenceLength + 1);
        var canonical = CanonicalSequence.Unmodified(longSequence);

        Assert.That(() => _serializer.Serialize(canonical), Throws.TypeOf<SequenceConversionException>());
    }

    [Test]
    public void Serialize_InvalidAminoAcid_ReturnsNullInReturnNullMode()
    {
        var canonical = CanonicalSequence.Unmodified("AX");
        var warnings = new ConversionWarnings();

        var result = _serializer.Serialize(canonical, warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(result, Is.Null);
        Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.InvalidSequence));
    }

    [Test]
    public void Serialize_UsesUnimodMappingsForResidueModifications()
    {
        var canonical = new CanonicalSequenceBuilder("M")
            .AddResidueModification(0, "Oxidation on M", unimodId: 35)
            .Build();

        var result = _serializer.Serialize(canonical);

        Assert.That(result, Is.EqualTo("-m_"));
    }

    [Test]
    public void Serialize_MassFallbackMappingHandlesRetentionMetadata()
    {
        var canonical = new CanonicalSequenceBuilder("AK")
            .WithSourceFormat("RetentionTimePrediction")
            .AddResidueModification(1, "RetentionAcetyl", mass: 42.01)
            .Build();

        var result = _serializer.Serialize(canonical);

        Assert.That(result, Is.EqualTo("-Aa_"));
    }

    [Test]
    public void Serialize_EmitsNterminusTokenForAcetylation()
    {
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddNTerminalModification("Common Biological:Acetylation on X", unimodId: 1)
            .Build();

        var result = _serializer.Serialize(canonical);

        Assert.That(result, Is.EqualTo("^PEPTIDE_"));
    }

    [Test]
    public void Serialize_UnrecognizedNTerm_ReturnNullAndRecordsWarning()
    {
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddNTerminalModification("RareNterm", mass: 13.3)
            .Build();
        var warnings = new ConversionWarnings();

        var result = _serializer.Serialize(canonical, warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(result, Is.Null);
        Assert.That(warnings.HasWarnings, Is.True);
        Assert.That(warnings.IncompatibleItems, Has.Some.EqualTo("N-terminal: RareNterm"));
    }

    [Test]
    public void Serialize_UnsupportedResidue_RemoveModeKeepsUnmodifiedSequence()
    {
        var canonical = CreateSequenceWithUnsupportedResidue();
        var warnings = new ConversionWarnings();

        var result = _serializer.Serialize(canonical, warnings, SequenceConversionHandlingMode.RemoveIncompatibleElements);

        var residueMod = canonical.Modifications.Single(m => m.PositionType == ModificationPositionType.Residue);

        Assert.That(result, Is.EqualTo("-PEPTIDE_"));
        Assert.That(warnings.HasWarnings, Is.True);
        Assert.That(warnings.IncompatibleItems, Has.Some.EqualTo(residueMod.ToString()));
    }

    [Test]
    public void Serialize_UnsupportedResidue_UsePrimarySequenceReturnsBase()
    {
        var canonical = CreateSequenceWithUnsupportedResidue();

        var result = _serializer.Serialize(canonical, null, SequenceConversionHandlingMode.UsePrimarySequence);

        Assert.That(result, Is.EqualTo("-PEPTIDE_"));
    }

    [Test]
    public void Serialize_UnsupportedResidue_ThrowModeFails()
    {
        var canonical = CreateSequenceWithUnsupportedResidue();

        Assert.That(() => _serializer.Serialize(canonical), Throws.TypeOf<SequenceConversionException>());
    }

    private static CanonicalSequence CreateSequenceWithUnsupportedResidue()
    {
        return new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(1, "UnsupportedRetention", mass: 123.45)
            .Build();
    }
}
