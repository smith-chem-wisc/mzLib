using NUnit.Framework;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion;

[TestFixture]
public class ChronologerSequenceFormatSchemaTests
{
    [Test]
    public void InstanceMetadataReflectsChronologerDefaults()
    {
        var schema = ChronologerSequenceFormatSchema.Instance;

        Assert.That(schema.FormatName, Is.EqualTo("Chronologer"));
        Assert.That(schema.ModOpenBracket, Is.EqualTo('['));
        Assert.That(schema.ModCloseBracket, Is.EqualTo(']'));
        Assert.That(schema.NTermSeparator, Is.Null);
        Assert.That(schema.CTermSeparator, Is.Null);
        Assert.That(schema.SupportsNTerminalMods, Is.False);
        Assert.That(schema.SupportsCTerminalMods, Is.False);
        Assert.That(schema.ToString(), Is.EqualTo("Chronologer"));
    }

    [Test]
    public void TerminusTokensMatchChronologerSpecification()
    {
        Assert.That(ChronologerSequenceFormatSchema.FreeNTerminus, Is.EqualTo('-'));
        Assert.That(ChronologerSequenceFormatSchema.CTerminus, Is.EqualTo('_'));
        Assert.That(ChronologerSequenceFormatSchema.NTermAcetyl, Is.EqualTo('^'));
        Assert.That(ChronologerSequenceFormatSchema.NTermPyroGlu, Is.EqualTo(')'));
        Assert.That(ChronologerSequenceFormatSchema.NTermCyclizedCamCys, Is.EqualTo('('));
        Assert.That(ChronologerSequenceFormatSchema.NTermTMT0, Is.EqualTo('&'));
        Assert.That(ChronologerSequenceFormatSchema.NTermTMT10, Is.EqualTo('*'));
    }

    [Test]
    public void AlphabetAndLengthConstraintsAreEnforced()
    {
        Assert.That(ChronologerSequenceFormatSchema.MaxSequenceLength, Is.EqualTo(50));
        Assert.That(ChronologerSequenceFormatSchema.EncodedLength, Is.EqualTo(52));
        Assert.That(ChronologerSequenceFormatSchema.CanonicalAminoAcids, Is.EqualTo("ACDEFGHIKLMNPQRSTVWY"));
        Assert.That(ChronologerSequenceFormatSchema.CanonicalAminoAcids.Length, Is.EqualTo(20));
    }
}
