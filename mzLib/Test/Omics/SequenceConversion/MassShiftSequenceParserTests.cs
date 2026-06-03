using NUnit.Framework;
using Omics.SequenceConversion;
using System.Globalization;
using System.Threading;

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
    public void CanParse_ReturnsFalseForWhitespaceInput()
    {
        var parser = new MassShiftSequenceParser();

        Assert.That(parser.CanParse("   "), Is.False);
    }

    [Test]
    public void CanParse_ReturnsFalseWhenClosingBracketOutnumbersOpening()
    {
        var parser = new MassShiftSequenceParser();

        Assert.That(parser.CanParse("[PEP]]"), Is.False);
    }

    [Test]
    public void CanParse_ReturnsFalseWhenNoNumericContent()
    {
        var parser = new MassShiftSequenceParser();

        Assert.That(parser.CanParse("PEP[NotMass]TIDE"), Is.False);
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

    [Test]
    public void Parse_InvalidMassShiftTokenWithThrowMode_Throws()
    {
        var parser = new MassShiftSequenceParser();

        Assert.That(
            () => parser.Parse("PEP[abc]TIDE", new ConversionWarnings(), SequenceConversionHandlingMode.ThrowException),
            Throws.TypeOf<SequenceConversionException>()
                .With.Property(nameof(SequenceConversionException.FailureReason))
                .EqualTo(ConversionFailureReason.UnknownFormat));
    }

    [Test]
    public void Parse_CTerminalSeparatorWithoutBracket_ReturnsNullAndFailure()
    {
        var parser = new MassShiftSequenceParser();
        var warnings = new ConversionWarnings();

        var result = parser.Parse("PEPTIDE-", warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(result, Is.Null);
        Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.UnknownFormat));
        Assert.That(warnings.Errors, Is.Not.Empty);
    }

    [Test]
    public void Parse_WithNTerminalResidueAndCTerminalMods_ReturnsCanonicalSequence()
    {
        var parser = new MassShiftSequenceParser();

        var sequence = parser.Parse("[+42.011]PEP[+15.995]TIDE-[+0.984]", null, SequenceConversionHandlingMode.ThrowException);

        Assert.That(sequence, Is.Not.Null);
        var canonical = sequence!.Value;

        Assert.That(canonical.BaseSequence, Is.EqualTo("PEPTIDE"));
        Assert.That(canonical.ModificationCount, Is.EqualTo(3));

        var nTerm = canonical.NTerminalModification;
        Assert.That(nTerm, Is.Not.Null);
        Assert.That(nTerm!.Value.PositionType, Is.EqualTo(ModificationPositionType.NTerminus));
        Assert.That(nTerm.Value.OriginalRepresentation, Is.EqualTo("+42.011"));

        var residueMod = canonical.GetModificationAt(2);
        Assert.That(residueMod, Is.Not.Null);
        Assert.That(residueMod!.Value.PositionType, Is.EqualTo(ModificationPositionType.Residue));
        Assert.That(residueMod.Value.TargetResidue, Is.EqualTo('P'));
        Assert.That(residueMod.Value.OriginalRepresentation, Is.EqualTo("+15.995"));

        var cTerm = canonical.CTerminalModification;
        Assert.That(cTerm, Is.Not.Null);
        Assert.That(cTerm!.Value.PositionType, Is.EqualTo(ModificationPositionType.CTerminus));
        Assert.That(cTerm.Value.TargetResidue, Is.EqualTo('E'));
        Assert.That(cTerm.Value.OriginalRepresentation, Is.EqualTo("+0.984"));
    }

    [Test]
    public void Parse_DoubleParsingFailureRespectsHandlingMode()
    {
        var parser = new MassShiftSequenceParser();
        var warnings = new ConversionWarnings();
        var originalCulture = CultureInfo.CurrentCulture;
        try
        {
            Thread.CurrentThread.CurrentCulture = new CultureInfo("fr-FR");

            var result = parser.Parse("PEP[+15.995]TIDE", warnings, SequenceConversionHandlingMode.ReturnNull);

            Assert.That(result, Is.Null);
            Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.UnknownFormat));
            Assert.That(warnings.HasErrors, Is.True);
        }
        finally
        {
            Thread.CurrentThread.CurrentCulture = originalCulture;
        }
    }
}
