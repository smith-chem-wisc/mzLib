using NUnit.Framework;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion;

[TestFixture]
public class SequenceParserBaseTests
{
    [Test]
    public void ParseMzLibFormatReturnsCanonicalSequenceWithTerminals()
    {
        var parser = new MzLibSequenceParser();
        var canonical = parser.Parse(
            "[Common Biological:Acetylation on S]SPEPTK[Common Fixed:Carbamidomethyl on K]IDE-[Common Artifact:Amidation on E]");

        Assert.That(canonical.HasValue, Is.True);
        var sequence = canonical.Value;

        Assert.That(sequence.BaseSequence, Is.EqualTo("SPEPTKIDE"));
        Assert.That(sequence.ModificationCount, Is.EqualTo(3));

        var nTerm = sequence.NTerminalModification;
        Assert.That(nTerm, Is.Not.Null);
        Assert.That(nTerm!.Value.TargetResidue, Is.EqualTo('S'));
        Assert.That(nTerm.Value.MzLibId, Is.EqualTo("Common Biological:Acetylation on S"));

        var residueMod = sequence.GetModificationAt(5);
        Assert.That(residueMod, Is.Not.Null);
        Assert.That(residueMod!.Value.TargetResidue, Is.EqualTo('K'));
        Assert.That(residueMod.Value.MzLibId, Is.EqualTo("Common Fixed:Carbamidomethyl on K"));

        var cTerm = sequence.CTerminalModification;
        Assert.That(cTerm, Is.Not.Null);
        Assert.That(cTerm!.Value.TargetResidue, Is.EqualTo('E'));
        Assert.That(cTerm.Value.MzLibId, Is.EqualTo("Common Artifact:Amidation on E"));
    }

    [Test]
    public void ParseRecordsWarningsForUnexpectedCharacters()
    {
        var parser = new TestSequenceParser(new TestSequenceFormatSchema("Test"));
        var warnings = new ConversionWarnings();

        var canonical = parser.Parse("A$B[mod]C", warnings);

        Assert.That(canonical.HasValue, Is.True);
        Assert.That(warnings.Warnings, Has.Count.EqualTo(1));
        Assert.That(warnings.Warnings[0], Does.Contain("Unexpected character '$'"));

        var sequence = canonical.Value;
        Assert.That(sequence.BaseSequence, Is.EqualTo("ABC"));
        var residueMod = sequence.GetModificationAt(1);
        Assert.That(residueMod, Is.Not.Null);
        Assert.That(residueMod!.Value.OriginalRepresentation, Is.EqualTo("mod"));
    }

    [Test]
    public void ParseEnforcesNTermSeparatorWhenSchemaRequiresOne()
    {
        var schema = new TestSequenceFormatSchema("SeparatorTest", nTermSeparator: "!");
        var parser = new TestSequenceParser(schema);
        var warnings = new ConversionWarnings();

        var result = parser.Parse("[mod]ABC", warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(result, Is.Null);
        Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.UnknownFormat));
        Assert.That(warnings.Errors, Has.Some.Contains("Expected N-terminal separator '!"));
    }

    private sealed class TestSequenceParser : SequenceParserBase
    {
        private readonly SequenceFormatSchema _schema;

        public TestSequenceParser(SequenceFormatSchema schema)
        {
            _schema = schema;
        }

        public override string FormatName => _schema.FormatName;

        public override SequenceFormatSchema Schema => _schema;

        public override bool CanParse(string input) => true;

        protected override CanonicalModification? ParseModificationString(
            string modString,
            ModificationPositionType positionType,
            int? residueIndex,
            char? targetResidue,
            ConversionWarnings warnings,
            SequenceConversionHandlingMode mode)
        {
            return new CanonicalModification(
                positionType,
                residueIndex,
                targetResidue,
                modString);
        }
    }

    private sealed class TestSequenceFormatSchema : SequenceFormatSchema
    {
        private readonly string _name;

        public TestSequenceFormatSchema(string name, string? nTermSeparator = "", string? cTermSeparator = "-")
            : base('[', ']', nTermSeparator, cTermSeparator)
        {
            _name = name;
        }

        public override string FormatName => _name;
    }
}
