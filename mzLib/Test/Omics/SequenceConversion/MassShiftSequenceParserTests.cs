using NUnit.Framework;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion;

[TestFixture]
public class MassShiftSequenceParserTests
{
    [Test]
    public void CanParse_ReturnsFalseWithoutBrackets()
    {
        var parser = new MassShiftSequenceParser();

        Assert.That(parser.CanParse("PEPTIDE"), Is.False);
    }

    [Test]
    public void CanParse_ReturnsFalseWhenBracketsUnbalanced()
    {
        var parser = new MassShiftSequenceParser();

        Assert.That(parser.CanParse("PEP[+15.99TIDE"), Is.False);
    }

    [Test]
    public void Parse_InvalidMassShiftToken_ReturnsNullAndFailure()
    {
        var parser = new MassShiftSequenceParser();
        var warnings = new ConversionWarnings();

        var result = parser.Parse("PEP[abc]TIDE", warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(result, Is.Null);
        Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.UnknownFormat));
    }

    [Test]
    public void Parse_ValidMassShift_ReturnsCanonicalSequence()
    {
        var parser = new MassShiftSequenceParser();

        var result = parser.Parse("PEP[+15.995]TIDE", null, SequenceConversionHandlingMode.ThrowException);

        Assert.That(result, Is.Not.Null);
        Assert.That(result!.Value.BaseSequence, Is.EqualTo("PEPTIDE"));
        Assert.That(result.Value.ModificationCount, Is.EqualTo(1));
    }
}
