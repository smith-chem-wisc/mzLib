using NUnit.Framework;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion;

[TestFixture]
public class UnimodSequenceParserTests
{
    [Test]
    public void CanParse_ReturnsTrueForLabeledUnimodTokens()
    {
        var parser = new UnimodSequenceParser();

        Assert.That(parser.CanParse("PEPC[UNIMOD:4]IDE"), Is.True);
    }

    [Test]
    public void CanParse_ReturnsFalseForMzLibStyleTokens()
    {
        var parser = new UnimodSequenceParser();

        Assert.That(parser.CanParse("PEPC[Common Fixed:Carbamidomethyl on C]IDE"), Is.False);
    }

    [Test]
    public void Parse_LabeledUnimodToken_SetsUnimodId()
    {
        var parser = new UnimodSequenceParser();

        var result = parser.Parse("PEPC[UNIMOD:4]IDE");

        Assert.That(result, Is.Not.Null);
        var canonical = result!.Value;
        var modification = canonical.GetModificationAt(3);
        Assert.That(modification, Is.Not.Null);
        Assert.That(modification!.Value.UnimodId, Is.EqualTo(4));
        Assert.That(modification.Value.TargetResidue, Is.EqualTo('C'));
    }

    [Test]
    public void Parse_UnlabeledUnimodToken_AcceptsIntegerId()
    {
        var parser = new UnimodSequenceParser();

        var result = parser.Parse("PEPC[4]IDE");

        Assert.That(result, Is.Not.Null);
        var canonical = result!.Value;
        var modification = canonical.GetModificationAt(3);
        Assert.That(modification, Is.Not.Null);
        Assert.That(modification!.Value.UnimodId, Is.EqualTo(4));
    }

    [Test]
    public void MzLibParser_CanParse_ReturnsFalseForLabeledUnimodTokens()
    {
        var parser = new MzLibSequenceParser();

        Assert.That(parser.CanParse("PEPC[UNIMOD:4]IDE"), Is.False);
    }
}
