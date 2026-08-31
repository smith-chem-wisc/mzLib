using Chemistry;
using NUnit.Framework;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion;

[TestFixture]
public class SequenceSerializerBaseTests
{
    private TestSequenceSerializer _serializer = null!;

    [SetUp]
    public void SetUp()
    {
        _serializer = new TestSequenceSerializer();
    }

    [Test]
    public void Serialize_IncludesTerminalAndResidueRepresentations()
    {
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddNTerminalModification("Nmod")
            .AddResidueModification(2, "Residue")
            .AddCTerminalModification("Cmod")
            .Build();

        var result = _serializer.Serialize(canonical);

        Assert.That(result, Is.EqualTo("[Nmod]-PEP[Residue]TIDE-[Cmod]"));
    }

    [Test]
    public void Serialize_SkipsUnsupportedResidueAndRecordsWarnings()
    {
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(3, TestSequenceSerializer.UnsupportedToken)
            .Build();
        var warnings = new ConversionWarnings();

        var result = _serializer.Serialize(canonical, warnings, SequenceConversionHandlingMode.RemoveIncompatibleElements);

        Assert.That(result, Is.EqualTo("PEPTIDE"));
        Assert.That(warnings.HasWarnings, Is.True);
        Assert.That(warnings.HasIncompatibleItems, Is.True);
        Assert.That(warnings.IncompatibleItems, Has.Some.EqualTo(TestSequenceSerializer.UnsupportedToken));
    }

    [Test]
    public void Serialize_ReturnsNullWhenUnsupportedInReturnNullMode()
    {
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(3, TestSequenceSerializer.UnsupportedToken)
            .Build();
        var warnings = new ConversionWarnings();

        var result = _serializer.Serialize(canonical, warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(result, Is.Null);
        Assert.That(warnings.HasWarnings, Is.True);
    }

    [Test]
    public void Serialize_EnrichesResolvableModifications()
    {
        // EnrichModificationsIfNeeded enriches mod metadata but preserves OriginalRepresentation
        // This test verifies that enrichment runs without error and still outputs the original token
        var lookup = new TestEnrichingLookup();
        var serializer = new TestSequenceSerializer(lookup);
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(1, TestSequenceSerializer.ResolvableToken)
            .Build();

        var result = serializer.Serialize(canonical);

        // After enrichment, the original representation is preserved in the output
        Assert.That(result, Does.Contain("[" + TestSequenceSerializer.ResolvableToken + "]"));
    }

    [Test]
    public void Serialize_EmptySequence_ReturnsNullAndRecordsFailure()
    {
        var warnings = new ConversionWarnings();

        var result = _serializer.Serialize(CanonicalSequence.Empty, warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(result, Is.Null);
        Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.InvalidSequence));
    }

    private sealed class TestSequenceSerializer : SequenceSerializerBase
    {
        public const string UnsupportedToken = "UnsupportedMod";
        public const string ResolvableToken = "ResolvableMod";

        private static readonly TestSequenceSchema SchemaInstance = new();

        public TestSequenceSerializer(IModificationLookup? lookup = null)
            : base(lookup)
        {
        }

        public override string FormatName => "TestFormat";
        public override SequenceFormatSchema Schema => SchemaInstance;
        public override bool CanSerialize(CanonicalSequence sequence) => true;
        public override bool ShouldResolveMod(CanonicalModification mod)
            => mod.OriginalRepresentation == ResolvableToken;

        protected override string? GetModificationString(CanonicalModification mod, ConversionWarnings warnings, SequenceConversionHandlingMode mode)
        {
            if (mod.OriginalRepresentation == UnsupportedToken)
            {
                warnings.AddWarning("Unsupported metadata encountered");
                warnings.AddIncompatibleItem(mod.OriginalRepresentation);
                return null;
            }

            return mod.OriginalRepresentation;
        }
    }

    private sealed class TestSequenceSchema : SequenceFormatSchema
    {
        public TestSequenceSchema()
            : base('[', ']', "-", "-")
        {
        }

        public override string FormatName => "TestSchema";
    }

    private sealed class TestEnrichingLookup : IModificationLookup
    {
        public string Name => "TestLookup";

        public CanonicalModification? TryResolve(CanonicalModification mod)
        {
            if (mod.OriginalRepresentation == TestSequenceSerializer.ResolvableToken)
            {
                return mod with { OriginalRepresentation = "resolved:" + mod.OriginalRepresentation };
            }

            return null;
        }

        public CanonicalModification? TryResolve(string originalRepresentation, char? targetResidue = null, ChemicalFormula? chemicalFormula = null, ModificationPositionType? positionType = null)
            => null;
    }
}
