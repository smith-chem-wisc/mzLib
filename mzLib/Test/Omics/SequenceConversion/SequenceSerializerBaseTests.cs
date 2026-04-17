using System;
using System.Collections.Generic;
using Chemistry;
using NUnit.Framework;
using Omics.Modifications;
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

    #region Serialize Tests

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
        var lookup = new TestEnrichingLookup();
        var serializer = new TestSequenceSerializer(lookup);
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(1, TestSequenceSerializer.ResolvableToken)
            .Build();

        var result = serializer.Serialize(canonical);

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

    [Test]
    public void Serialize_UnexpectedException_ReturnsNullInReturnNullMode()
    {
        var serializer = new ControlledThrowSerializer(throwConversionException: false);
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(1, "AnyMod")
            .Build();
        var warnings = new ConversionWarnings();

        var result = serializer.Serialize(canonical, warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(result, Is.Null);
        Assert.That(warnings.HasFatalError, Is.True);
        Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.UnknownFormat));
    }

    [Test]
    public void Serialize_UnexpectedException_SetsFailureInRemoveIncompatibleMode()
    {
        var serializer = new ControlledThrowSerializer(throwConversionException: false);
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(1, "AnyMod")
            .Build();
        var warnings = new ConversionWarnings();

        var result = serializer.Serialize(canonical, warnings, SequenceConversionHandlingMode.RemoveIncompatibleElements);

        Assert.That(result, Is.Null);
        Assert.That(warnings.HasFatalError, Is.True);
    }

    [Test]
    public void Serialize_SequenceConversionException_Rethrown()
    {
        var serializer = new ControlledThrowSerializer(throwConversionException: true);
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(1, "AnyMod")
            .Build();

        var ex = Assert.Throws<SequenceConversionException>(() =>
            serializer.Serialize(canonical, new ConversionWarnings(), SequenceConversionHandlingMode.ThrowException));

        Assert.That(ex!.FailureReason, Is.EqualTo(ConversionFailureReason.IncompatibleModifications));
    }

    [Test]
    public void Serialize_NTermNullModString_SkipsInRemoveIncompatibleMode()
    {
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddNTerminalModification(TestSequenceSerializer.UnsupportedToken)
            .Build();
        var warnings = new ConversionWarnings();

        var result = _serializer.Serialize(canonical, warnings, SequenceConversionHandlingMode.RemoveIncompatibleElements);

        Assert.That(result, Is.EqualTo("PEPTIDE"));
        Assert.That(warnings.HasWarnings, Is.True);
    }

    [Test]
    public void Serialize_NTermNullModString_ReturnsNullInReturnNullMode()
    {
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddNTerminalModification(TestSequenceSerializer.UnsupportedToken)
            .Build();
        var warnings = new ConversionWarnings();

        var result = _serializer.Serialize(canonical, warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(result, Is.Null);
    }

    [Test]
    public void Serialize_CTermNullModString_SkipsInRemoveIncompatibleMode()
    {
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddCTerminalModification(TestSequenceSerializer.UnsupportedToken)
            .Build();
        var warnings = new ConversionWarnings();

        var result = _serializer.Serialize(canonical, warnings, SequenceConversionHandlingMode.RemoveIncompatibleElements);

        Assert.That(result, Is.EqualTo("PEPTIDE"));
        Assert.That(warnings.HasWarnings, Is.True);
    }

    [Test]
    public void Serialize_CTermNullModString_ReturnsNullInReturnNullMode()
    {
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddCTerminalModification(TestSequenceSerializer.UnsupportedToken)
            .Build();
        var warnings = new ConversionWarnings();

        var result = _serializer.Serialize(canonical, warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(result, Is.Null);
    }

    #endregion

    #region EnrichModificationsIfNeeded Tests (via Serialize)

    [Test]
    public void EnrichModificationsIfNeeded_NoLookup_OutputUnchanged()
    {
        var serializer = new TestSequenceSerializer(lookup: null);
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(1, TestSequenceSerializer.ResolvableToken)
            .Build();

        var result = serializer.Serialize(canonical);

        Assert.That(result, Does.Contain("[" + TestSequenceSerializer.ResolvableToken + "]"));
    }

    [Test]
    public void EnrichModificationsIfNeeded_ShouldResolveFalse_LookupNotCalled()
    {
        var recordingLookup = new RecordingLookup();
        var serializer = new TestSequenceSerializer(recordingLookup);
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(1, "NonResolvableToken")
            .Build();

        var result = serializer.Serialize(canonical);

        Assert.That(result, Does.Contain("[NonResolvableToken]"));
        Assert.That(recordingLookup.TryResolveCallCount, Is.EqualTo(0));
    }

    [Test]
    public void EnrichModificationsIfNeeded_LookupFails_OutputUnchanged()
    {
        var recordingLookup = new RecordingLookup();
        var serializer = new TestSequenceSerializer(recordingLookup);
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(1, TestSequenceSerializer.ResolvableToken)
            .Build();

        var result = serializer.Serialize(canonical);

        Assert.That(result, Does.Contain("[" + TestSequenceSerializer.ResolvableToken + "]"));
        Assert.That(recordingLookup.TryResolveCallCount, Is.EqualTo(1));
    }

    #endregion

    #region ToOneIsNterminusModificationDictionary Tests

    [Test]
    public void ToOneIsNterminusModificationDictionary_EmptyBaseSequence_ThrowMode_ThrowsException()
    {
        var sequence = CanonicalSequence.Unmodified("");

        Assert.Throws<SequenceConversionException>(() =>
            _serializer.ToOneIsNterminusModificationDictionary(sequence, mode: SequenceConversionHandlingMode.ThrowException));
    }

    [Test]
    public void ToOneIsNterminusModificationDictionary_EmptyBaseSequence_NonThrowMode_ReturnsEmptyDict()
    {
        var sequence = CanonicalSequence.Unmodified("");
        var warnings = new ConversionWarnings();

        var result = _serializer.ToOneIsNterminusModificationDictionary(
            sequence, warnings: warnings, mode: SequenceConversionHandlingMode.RemoveIncompatibleElements);

        Assert.That(result, Is.Empty);
        Assert.That(warnings.HasFatalError, Is.True);
        Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.InvalidSequence));
    }

    [Test]
    public void ToOneIsNterminusModificationDictionary_UnmodifiedSequence_ReturnsEmptyDict()
    {
        var sequence = CanonicalSequence.Unmodified("PEPTIDE");

        var result = _serializer.ToOneIsNterminusModificationDictionary(sequence);

        Assert.That(result, Is.Empty);
    }

    [Test]
    public void ToOneIsNterminusModificationDictionary_ResolvedResidueMod_MapsToCorrectIndex()
    {
        var oxMod = Mods.AllKnownProteinModsDictionary["Oxidation on M"];
        var sequence = new CanonicalSequenceBuilder("PEPTMIDE")
            .AddResidueModification(4, "Oxidation", mzLibModification: oxMod)
            .Build();

        var result = _serializer.ToOneIsNterminusModificationDictionary(sequence);

        Assert.That(result, Contains.Key(6));
        Assert.That(result[6], Is.SameAs(oxMod));
    }

    [Test]
    public void ToOneIsNterminusModificationDictionary_ResolvedNTerminalMod_MapsToOne()
    {
        var acetylMod = GetNTerminalModification();
        var sequence = new CanonicalSequenceBuilder("PEPTIDE")
            .AddNTerminalModification("Acetyl", mzLibModification: acetylMod)
            .Build();

        var result = _serializer.ToOneIsNterminusModificationDictionary(sequence);

        Assert.That(result, Contains.Key(1));
        Assert.That(result[1], Is.SameAs(acetylMod));
    }

    [Test]
    public void ToOneIsNterminusModificationDictionary_ResolvedCTerminalMod_MapsToLengthPlusTwo()
    {
        var cTermMod = GetCTerminalModification();
        var sequence = new CanonicalSequenceBuilder("PEPTIDE")
            .AddCTerminalModification("CTermMod", mzLibModification: cTermMod)
            .Build();

        var result = _serializer.ToOneIsNterminusModificationDictionary(sequence);

        Assert.That(result, Contains.Key(9));
        Assert.That(result[9], Is.SameAs(cTermMod));
    }

    [Test]
    public void ToOneIsNterminusModificationDictionary_MultipleResolvedMods_AllMappedCorrectly()
    {
        var nTermMod = GetNTerminalModification();
        var oxMod = Mods.AllKnownProteinModsDictionary["Oxidation on M"];
        var cTermMod = GetCTerminalModification();

        var sequence = new CanonicalSequenceBuilder("PEPTMIDE")
            .AddNTerminalModification("Acetyl", mzLibModification: nTermMod)
            .AddResidueModification(4, "Oxidation", mzLibModification: oxMod)
            .AddCTerminalModification("CTerm", mzLibModification: cTermMod)
            .Build();

        var result = _serializer.ToOneIsNterminusModificationDictionary(sequence);

        Assert.That(result, Has.Count.EqualTo(3));
        Assert.That(result, Contains.Key(1));
        Assert.That(result[1], Is.SameAs(nTermMod));
        Assert.That(result, Contains.Key(6));
        Assert.That(result[6], Is.SameAs(oxMod));
        Assert.That(result, Contains.Key(10));
        Assert.That(result[10], Is.SameAs(cTermMod));
    }

    [Test]
    public void ToOneIsNterminusModificationDictionary_UnresolvedMod_ThrowMode_ThrowsException()
    {
        var sequence = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(1, "UnknownMod")
            .Build();

        var ex = Assert.Throws<SequenceConversionException>(() =>
            _serializer.ToOneIsNterminusModificationDictionary(
                sequence, mode: SequenceConversionHandlingMode.ThrowException));

        Assert.That(ex!.FailureReason, Is.EqualTo(ConversionFailureReason.IncompatibleModifications));
    }

    [Test]
    public void ToOneIsNterminusModificationDictionary_UnresolvedMod_WarnMode_SkipsAndWarns()
    {
        var oxMod = Mods.AllKnownProteinModsDictionary["Oxidation on M"];
        var sequence = new CanonicalSequenceBuilder("PEPTMIDE")
            .AddResidueModification(1, "UnknownMod")
            .AddResidueModification(4, "Oxidation", mzLibModification: oxMod)
            .Build();
        var warnings = new ConversionWarnings();

        var result = _serializer.ToOneIsNterminusModificationDictionary(
            sequence, warnings: warnings, mode: SequenceConversionHandlingMode.RemoveIncompatibleElements);

        Assert.That(result, Has.Count.EqualTo(1));
        Assert.That(result, Contains.Key(6));
        Assert.That(warnings.HasWarnings, Is.True);
        Assert.That(warnings.HasIncompatibleItems, Is.True);
    }

    [Test]
    public void ToOneIsNterminusModificationDictionary_DuplicatePosition_ThrowMode_ThrowsException()
    {
        var oxMod = Mods.AllKnownProteinModsDictionary["Oxidation on M"];
        var carbMod = Mods.AllKnownProteinModsDictionary["Carbamidomethyl on C"];
        var mod1 = CanonicalModification.AtResidue(2, 'P', "Mod1", mzLibModification: oxMod);
        var mod2 = CanonicalModification.AtResidue(2, 'P', "Mod2", mzLibModification: carbMod);
        var sequence = new CanonicalSequenceBuilder("PEPTIDE")
            .AddModification(mod1)
            .AddModification(mod2)
            .Build();

        var ex = Assert.Throws<SequenceConversionException>(() =>
            _serializer.ToOneIsNterminusModificationDictionary(
                sequence, mode: SequenceConversionHandlingMode.ThrowException));

        Assert.That(ex!.FailureReason, Is.EqualTo(ConversionFailureReason.InvalidSequence));
    }

    [Test]
    public void ToOneIsNterminusModificationDictionary_DuplicatePosition_WarnMode_KeepsFirst()
    {
        var oxMod = Mods.AllKnownProteinModsDictionary["Oxidation on M"];
        var carbMod = Mods.AllKnownProteinModsDictionary["Carbamidomethyl on C"];
        var mod1 = CanonicalModification.AtResidue(2, 'P', "Mod1", mzLibModification: oxMod);
        var mod2 = CanonicalModification.AtResidue(2, 'P', "Mod2", mzLibModification: carbMod);
        var sequence = new CanonicalSequenceBuilder("PEPTIDE")
            .AddModification(mod1)
            .AddModification(mod2)
            .Build();
        var warnings = new ConversionWarnings();

        var result = _serializer.ToOneIsNterminusModificationDictionary(
            sequence, warnings: warnings, mode: SequenceConversionHandlingMode.RemoveIncompatibleElements);

        Assert.That(result, Has.Count.EqualTo(1));
        Assert.That(result[4], Is.SameAs(oxMod));
        Assert.That(warnings.HasWarnings, Is.True);
    }

    [Test]
    public void ToOneIsNterminusModificationDictionary_InvalidPositionType_ThrowsException()
    {
        var invalidMod = new CanonicalModification((ModificationPositionType)999, 0, 'P', "InvalidPosType");
        var sequence = new CanonicalSequenceBuilder("PEPTIDE")
            .AddModification(invalidMod)
            .Build();

        var ex = Assert.Throws<SequenceConversionException>(() =>
            _serializer.ToOneIsNterminusModificationDictionary(sequence));

        Assert.That(ex!.Message, Does.Contain("Unsupported modification position type"));
    }

    [Test]
    public void ToOneIsNterminusModificationDictionary_ResidueMod_NullIndex_ThrowsException()
    {
        var nullIndexMod = new CanonicalModification(ModificationPositionType.Residue, null, 'P', "NullIndex");
        var sequence = new CanonicalSequenceBuilder("PEPTIDE")
            .AddModification(nullIndexMod)
            .Build();

        var ex = Assert.Throws<SequenceConversionException>(() =>
            _serializer.ToOneIsNterminusModificationDictionary(sequence));

        Assert.That(ex!.Message, Does.Contain("invalid"));
    }

    [Test]
    public void ToOneIsNterminusModificationDictionary_ResidueMod_NegativeIndex_ThrowsException()
    {
        var negIndexMod = new CanonicalModification(ModificationPositionType.Residue, -1, 'P', "NegIndex");
        var sequence = new CanonicalSequenceBuilder("PEPTIDE")
            .AddModification(negIndexMod)
            .Build();

        var ex = Assert.Throws<SequenceConversionException>(() =>
            _serializer.ToOneIsNterminusModificationDictionary(sequence));

        Assert.That(ex!.Message, Does.Contain("invalid"));
    }

    [Test]
    public void ToOneIsNterminusModificationDictionary_ResidueMod_IndexOutOfRange_ThrowsException()
    {
        var outOfRangeMod = new CanonicalModification(ModificationPositionType.Residue, 100, 'P', "OutOfRange");
        var sequence = new CanonicalSequenceBuilder("PEPTIDE")
            .AddModification(outOfRangeMod)
            .Build();

        var ex = Assert.Throws<SequenceConversionException>(() =>
            _serializer.ToOneIsNterminusModificationDictionary(sequence));

        Assert.That(ex!.Message, Does.Contain("invalid"));
    }

    [Test]
    public void ToOneIsNterminusModificationDictionary_UnexpectedException_NonThrowMode_ReturnsEmptyDict()
    {
        var throwingLookup = new ThrowingLookup();
        var serializer = new TestSequenceSerializer(throwingLookup);
        var sequence = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(1, "UnresolvedMod")
            .Build();
        var warnings = new ConversionWarnings();

        var result = serializer.ToOneIsNterminusModificationDictionary(
            sequence, warnings: warnings, mode: SequenceConversionHandlingMode.RemoveIncompatibleElements);

        Assert.That(result, Is.Empty);
        Assert.That(warnings.HasFatalError, Is.True);
        Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.UnknownFormat));
    }

    #endregion

    #region ResolveModificationsForProjection Tests (via ToOneIsNterminusModificationDictionary)

    [Test]
    public void ResolveModificationsForProjection_NoLookups_ModStaysUnresolved_ThrowsInThrowMode()
    {
        var serializer = new TestSequenceSerializer(lookup: null);
        var sequence = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(1, "UnresolvedMod")
            .Build();

        Assert.Throws<SequenceConversionException>(() =>
            serializer.ToOneIsNterminusModificationDictionary(sequence, mode: SequenceConversionHandlingMode.ThrowException));
    }

    [Test]
    public void ResolveModificationsForProjection_AlreadyResolved_SkipsLookup()
    {
        var throwingLookup = new ThrowingLookup();
        var serializer = new TestSequenceSerializer(throwingLookup);
        var oxMod = Mods.AllKnownProteinModsDictionary["Oxidation on M"];
        var sequence = new CanonicalSequenceBuilder("PEPTMIDE")
            .AddResidueModification(4, "Oxidation", mzLibModification: oxMod)
            .Build();

        var result = serializer.ToOneIsNterminusModificationDictionary(sequence);

        Assert.That(result, Contains.Key(6));
        Assert.That(result[6], Is.SameAs(oxMod));
    }

    [Test]
    public void ResolveModificationsForProjection_FallbackLookupResolves_UpdatesMod()
    {
        var serializer = new TestSequenceSerializer(lookup: null);
        var oxMod = Mods.AllKnownProteinModsDictionary["Oxidation on M"];
        var knownMods = new Dictionary<string, Modification> { { "Oxidation on M", oxMod } };
        var sequence = new CanonicalSequenceBuilder("PEPTMIDE")
            .AddResidueModification(4, "Oxidation on M", mzLibId: "Oxidation on M")
            .Build();

        var result = serializer.ToOneIsNterminusModificationDictionary(sequence, knownMods: knownMods);

        Assert.That(result, Contains.Key(6));
    }

    [Test]
    public void ResolveModificationsForProjection_PrimaryLookupResolves_WhenFallbackFails()
    {
        var oxMod = Mods.AllKnownProteinModsDictionary["Oxidation on M"];
        var selectiveLookup = new SelectiveResolveLookup(
            new Dictionary<string, Modification> { { "Oxidation", oxMod } });
        var serializer = new TestSequenceSerializer(selectiveLookup);
        var emptyKnownMods = new Dictionary<string, Modification>();
        var sequence = new CanonicalSequenceBuilder("PEPTMIDE")
            .AddResidueModification(4, "Oxidation")
            .Build();

        var result = serializer.ToOneIsNterminusModificationDictionary(sequence, knownMods: emptyKnownMods);

        Assert.That(result, Contains.Key(6));
    }

    [Test]
    public void ResolveModificationsForProjection_BothLookupsFail_ModStaysUnresolved()
    {
        var serializer = new TestSequenceSerializer(lookup: null);
        var emptyKnownMods = new Dictionary<string, Modification>();
        var sequence = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(1, "CompletelyUnknown")
            .Build();
        var warnings = new ConversionWarnings();

        var result = serializer.ToOneIsNterminusModificationDictionary(
            sequence, knownMods: emptyKnownMods, warnings: warnings,
            mode: SequenceConversionHandlingMode.RemoveIncompatibleElements);

        Assert.That(result, Is.Empty);
        Assert.That(warnings.HasWarnings, Is.True);
    }

    [Test]
    public void ResolveModificationsForProjection_FallbackResolvesBeforePrimaryLookup()
    {
        var oxMod = Mods.AllKnownProteinModsDictionary["Oxidation on M"];
        var carbMod = Mods.AllKnownProteinModsDictionary["Carbamidomethyl on C"];
        var primaryLookup = new SelectiveResolveLookup(
            new Dictionary<string, Modification> { { "Oxidation", carbMod } });
        var serializer = new TestSequenceSerializer(primaryLookup);
        var fallbackMods = new Dictionary<string, Modification> { { "Oxidation on M", oxMod } };
        var sequence = new CanonicalSequenceBuilder("PEPTMIDE")
            .AddResidueModification(4, "Oxidation on M", mzLibId: "Oxidation on M")
            .Build();

        var result = serializer.ToOneIsNterminusModificationDictionary(sequence, knownMods: fallbackMods);

        Assert.That(result, Contains.Key(6));
        Assert.That(result[6], Is.SameAs(oxMod));
    }

    #endregion

    #region Test Stubs

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

    private sealed class SelectiveResolveLookup : IModificationLookup
    {
        private readonly Dictionary<string, Modification> _resolved;

        public string Name => "Selective";

        public SelectiveResolveLookup(Dictionary<string, Modification> resolvedByOriginalRep)
        {
            _resolved = resolvedByOriginalRep;
        }

        public CanonicalModification? TryResolve(CanonicalModification mod)
        {
            if (_resolved.TryGetValue(mod.OriginalRepresentation, out var modification))
            {
                return mod.WithResolvedModification(modification, mod.ResidueIndex, mod.PositionType);
            }

            return null;
        }

        public CanonicalModification? TryResolve(string originalRepresentation, char? targetResidue = null,
            ChemicalFormula? chemicalFormula = null, ModificationPositionType? positionType = null)
        {
            if (_resolved.TryGetValue(originalRepresentation, out var modification))
            {
                return CanonicalModification.AtResidue(0, targetResidue ?? 'X', originalRepresentation,
                    mzLibModification: modification);
            }

            return null;
        }
    }

    private sealed class ThrowingLookup : IModificationLookup
    {
        public string Name => "Throwing";

        public CanonicalModification? TryResolve(CanonicalModification mod)
        {
            throw new InvalidOperationException("Lookup intentionally failed");
        }

        public CanonicalModification? TryResolve(string originalRepresentation, char? targetResidue = null,
            ChemicalFormula? chemicalFormula = null, ModificationPositionType? positionType = null)
        {
            throw new InvalidOperationException("Lookup intentionally failed");
        }
    }

    private sealed class RecordingLookup : IModificationLookup
    {
        public int TryResolveCallCount { get; private set; }
        public string Name => "Recording";

        public CanonicalModification? TryResolve(CanonicalModification mod)
        {
            TryResolveCallCount++;
            return null;
        }

        public CanonicalModification? TryResolve(string originalRepresentation, char? targetResidue = null,
            ChemicalFormula? chemicalFormula = null, ModificationPositionType? positionType = null)
        {
            TryResolveCallCount++;
            return null;
        }
    }

    private sealed class ControlledThrowSerializer : SequenceSerializerBase
    {
        private readonly bool _throwConversionException;
        private static readonly TestSequenceSchema SchemaInstance = new();

        public ControlledThrowSerializer(bool throwConversionException = false)
            : base(null)
        {
            _throwConversionException = throwConversionException;
        }

        public override string FormatName => "ControlledThrow";
        public override SequenceFormatSchema Schema => SchemaInstance;
        public override bool CanSerialize(CanonicalSequence sequence) => true;
        public override bool ShouldResolveMod(CanonicalModification mod) => false;

        protected override string? GetModificationString(CanonicalModification mod, ConversionWarnings warnings, SequenceConversionHandlingMode mode)
        {
            if (_throwConversionException)
            {
                throw new SequenceConversionException("test error", ConversionFailureReason.IncompatibleModifications);
            }

            throw new InvalidOperationException("unexpected error");
        }
    }

    #endregion

    #region Helpers

    private static Modification GetNTerminalModification()
    {
        ModificationMotif.TryGetMotif("X", out var motif);
        return new Modification("Acetyl", _modificationType: "Common Variable",
            _target: motif, _locationRestriction: "Peptide N-terminal.",
            _chemicalFormula: new ChemicalFormula(ChemicalFormula.ParseFormula("C2H2O1")),
            _monoisotopicMass: 42.010565);
    }

    private static Modification GetCTerminalModification()
    {
        ModificationMotif.TryGetMotif("X", out var motif);
        return new Modification("Amidated", _modificationType: "Common Variable",
            _target: motif, _locationRestriction: "Peptide C-terminal.",
            _chemicalFormula: new ChemicalFormula(ChemicalFormula.ParseFormula("H1N1")),
            _monoisotopicMass: -0.984016);
    }

    #endregion
}
