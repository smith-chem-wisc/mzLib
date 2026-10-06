using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Reflection;
using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using Readers.ProForma;

namespace Test.KoinaTests;

[TestFixture]
public class KoinaModelBaseTests
{
    private sealed class FakeSequenceConverter : ISequenceConverter
    {
        private readonly Func<string, CanonicalSequence?> _parse;
        private readonly Func<CanonicalSequence, string?> _serialize;

        public FakeSequenceConverter(Func<string, CanonicalSequence?> parse, Func<CanonicalSequence, string?> serialize)
        {
            _parse = parse;
            _serialize = serialize;
            Parser = new FakeSequenceParser(parse);
        }

        public string FormatName => "fake-fake";
        public string SourceFormatName => "fake";
        public string TargetFormatName => "fake";
        public ISequenceParser Parser { get; }
        public ISequenceSerializer Serializer => null!;

        public CanonicalSequence? Parse(string input, ConversionWarnings? warnings = null, SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
            => _parse(input);

        public string? Serialize(CanonicalSequence sequence, ConversionWarnings? warnings = null, SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
            => _serialize(sequence);

        public string? Convert(string input, ConversionWarnings? warnings = null, SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
        {
            var canonical = Parse(input, warnings, mode);
            return canonical.HasValue ? Serialize(canonical.Value, warnings, mode) : null;
        }
    }

    private sealed class BraceSchema() : SequenceFormatSchema('{', '}')
    {
        public override string FormatName => "braces";
    }

    private sealed class FakeSequenceParser : ISequenceParser
    {
        private readonly Func<string, CanonicalSequence?> _parse;

        public FakeSequenceParser(Func<string, CanonicalSequence?> parse, SequenceFormatSchema? schema = null)
        {
            _parse = parse;
            Schema = schema ?? MzLibSequenceFormatSchema.Instance;
        }

        public string FormatName => "fake";
        public SequenceFormatSchema Schema { get; }
        public bool CanParse(string input) => true;

        public CanonicalSequence? Parse(string input, ConversionWarnings? warnings = null, SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
            => _parse(input);
    }

    private sealed class KoinaModelHarness : KoinaModelBase<string, string>
    {
        public KoinaModelHarness(
            ISequenceConverter converter,
            SequenceConversionHandlingMode modHandlingMode = SequenceConversionHandlingMode.ReturnNull,
            IReadOnlySet<int>? allowedUnimodIds = null,
            bool acceptsAllUnimodModifications = false,
            IReadOnlySet<int>? requiredNTerminalUnimodIds = null)
            : base(converter)
        {
            ModHandlingMode = modHandlingMode;
            AllowedUnimodIds = allowedUnimodIds ?? new HashSet<int>();
            AcceptsAllUnimodModifications = acceptsAllUnimodModifications;
            RequiredNTerminalUnimodIds = requiredNTerminalUnimodIds;
        }

        public override string ModelName => "Harness";
        public override int MaxBatchSize => 10;
        public override int MaxNumberOfBatchesPerRequest { get; init; } = 1;
        public override int ThrottlingDelayInMilliseconds { get; init; } = 0;
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 0;
        public override SequenceConversionHandlingMode ModHandlingMode { get; init; }
        public override int MaxPeptideLength => 50;
        public override int MinPeptideLength => 1;
        public override IReadOnlySet<int> AllowedUnimodIds { get; }
        public override bool AcceptsAllUnimodModifications { get; }
        public override IReadOnlySet<int>? RequiredNTerminalUnimodIds { get; }

        protected override List<Dictionary<string, object>> ToBatchedRequests(List<string> validInputs)
        {
            return new List<Dictionary<string, object>>();
        }

        public string? TryClean(string sequence, out string? apiSequence, out WarningException? warning)
        {
            return TryCleanSequence(sequence, null, out apiSequence, out warning);
        }

        public string? TryCleanWithParser(string sequence, ISequenceParser? sourceParser, out string? apiSequence, out WarningException? warning)
        {
            return TryCleanSequence(sequence, sourceParser, out apiSequence, out warning);
        }

        public static ISequenceConverter BuildConverter(IReadOnlySet<int> allowedUnimodIds)
        {
            return CreateUnimodConverter(UnimodSequenceFormatSchema.Instance, allowedUnimodIds);
        }

        public static ISequenceConverter BuildAcceptAllConverter()
        {
            return CreateUnimodConverterAcceptAll(UnimodSequenceFormatSchema.Instance);
        }
    }

    [Test]
    public void Constructor_WithNullConverter_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => new KoinaModelHarness(null!));
    }

    [Test]
    public void TryCleanSequence_InvalidBaseSequence_ReturnsNull()
    {
        var model = new KoinaModelHarness(KoinaModelHarness.BuildConverter(new HashSet<int>()));

        var result = model.TryClean("PEP*TIDE", out var apiSequence, out var warning);

        Assert.That(result, Is.Null);
        Assert.That(apiSequence, Is.Null);
        Assert.That(warning?.Message, Does.Contain("Invalid base sequence 'PEP*TIDE'"));
    }

    // Everything outside the source format's complete modification brackets must be a residue the model allows, so
    // whitespace, an unknown character or an unbalanced bracket fails as a residue whichever parser reads it.
    [TestCase("PEPTIDE K", false)]
    [TestCase("PEP*TIDE", true)]
    [TestCase("PEPM[UNIMOD:35IDE", true)]
    public void TryCleanSequence_NonResidueOutsideModifications_IsRejected(string sequence, bool proForma)
    {
        var model = new KoinaModelHarness(KoinaModelHarness.BuildConverter(new HashSet<int> { 35 }), allowedUnimodIds: new HashSet<int> { 35 });

        var result = model.TryCleanWithParser(sequence, proForma ? ProFormaSequenceParser.Instance : null, out var apiSequence, out var warning);

        Assert.That(result, Is.Null);
        Assert.That(apiSequence, Is.Null);
        Assert.That(warning?.Message, Does.Contain("Invalid base sequence"));
    }

    [TestCase("[UNIMOD:737]-PEPTIDEK", true, "[UNIMOD:737]PEPTIDEK")]
    [TestCase("PEPTIDEK-[UNIMOD:2]", true, "PEPTIDEK-[UNIMOD:2]")]
    [TestCase("[Multiplex Label:TMT6-plex on X]PEPTIDEK", false, "[UNIMOD:737]PEPTIDEK")]
    [TestCase("PEPTIDEK-[Unimod:Amidated on X]", false, "PEPTIDEK-[UNIMOD:2]")]
    public void TryCleanSequence_TerminalModificationSeparators_AreNotResidues(string sequence, bool proForma, string expected)
    {
        var allowed = new HashSet<int> { 737, 2 };
        var model = new KoinaModelHarness(KoinaModelHarness.BuildConverter(allowed), allowedUnimodIds: allowed);

        var result = model.TryCleanWithParser(sequence, proForma ? ProFormaSequenceParser.Instance : null, out var apiSequence, out var warning);

        Assert.That(result, Is.Not.Null, warning?.Message);
        Assert.That(apiSequence, Is.EqualTo(expected));
    }

    // "none" = null (nothing required), "any" = empty (some N-terminal mod required), otherwise the required ids.
    private static IReadOnlySet<int>? RequiredIds(string required) => required switch
    {
        "none" => null,
        "any" => new HashSet<int>(),
        _ => required.Split(',').Select(int.Parse).ToHashSet()
    };

    [TestCase("none", "PEPTIDEK")]
    [TestCase("any", "[UNIMOD:214]-PEPTIDEK")]
    [TestCase("737", "[UNIMOD:737]-PEPTIDEK")]
    public void TryCleanSequence_RequiredNTerminalModification_AcceptsWhatItsSentinelAllows(string required, string sequence)
    {
        var allowed = new HashSet<int> { 737, 214 };
        var model = new KoinaModelHarness(KoinaModelHarness.BuildConverter(allowed), allowedUnimodIds: allowed,
            requiredNTerminalUnimodIds: RequiredIds(required));

        var result = model.TryCleanWithParser(sequence, ProFormaSequenceParser.Instance, out _, out var warning);

        Assert.That(result, Is.Not.Null, warning?.Message);
    }

    [TestCase("any", "PEPTIDEK")]
    [TestCase("737", "PEPTIDEK")]
    [TestCase("737", "[UNIMOD:214]-PEPTIDEK")]
    public void TryCleanSequence_RequiredNTerminalModification_RejectsWhatItsSentinelDoesNot(string required, string sequence)
    {
        var allowed = new HashSet<int> { 737, 214 };
        var model = new KoinaModelHarness(KoinaModelHarness.BuildConverter(allowed), allowedUnimodIds: allowed,
            requiredNTerminalUnimodIds: RequiredIds(required));

        var result = model.TryCleanWithParser(sequence, ProFormaSequenceParser.Instance, out var apiSequence, out var warning);

        Assert.That(result, Is.Null);
        Assert.That(apiSequence, Is.Null);
        Assert.That(warning?.Message, Does.Contain("N-terminal"));
    }

    [Test]
    public void TryCleanSequence_MissingRequiredNTerminalModification_FailsClosedInEveryMode(
        [Values(SequenceConversionHandlingMode.ReturnNull, SequenceConversionHandlingMode.RemoveIncompatibleElements,
            SequenceConversionHandlingMode.UsePrimarySequence)] SequenceConversionHandlingMode mode)
    {
        var allowed = new HashSet<int> { 737 };
        var model = new KoinaModelHarness(KoinaModelHarness.BuildConverter(allowed), mode, allowed, requiredNTerminalUnimodIds: allowed);

        var result = model.TryCleanWithParser("PEPTIDEK", ProFormaSequenceParser.Instance, out var apiSequence, out var warning);

        Assert.That(result, Is.Null);
        Assert.That(apiSequence, Is.Null);
        Assert.That(warning?.Message, Does.Contain("N-terminal"));
    }

    [Test]
    public void TryCleanSequence_MissingRequiredNTerminalModificationThrowMode_Throws()
    {
        var allowed = new HashSet<int> { 737 };
        var model = new KoinaModelHarness(KoinaModelHarness.BuildConverter(allowed),
            SequenceConversionHandlingMode.ThrowException, allowed, requiredNTerminalUnimodIds: allowed);

        Assert.That(() => model.TryCleanWithParser("PEPTIDEK", ProFormaSequenceParser.Instance, out _, out _),
            Throws.ArgumentException.With.Message.Contains("N-terminal"));
    }

    [Test]
    public void TryCleanSequence_RemovedRequiredNTerminalLabel_IsRejectedWithTheRemovalReason()
    {
        // UNIMOD:739 is outside the allow-list, so RemoveIncompatibleElements drops it, and the peptide is then
        // missing its required label. The warning must say both.
        var allowed = new HashSet<int> { 737 };
        var model = new KoinaModelHarness(KoinaModelHarness.BuildConverter(allowed),
            SequenceConversionHandlingMode.RemoveIncompatibleElements, allowed, requiredNTerminalUnimodIds: allowed);

        var result = model.TryCleanWithParser("[UNIMOD:739]-PEPTIDEK", ProFormaSequenceParser.Instance, out _, out var warning);

        Assert.That(result, Is.Null);
        Assert.That(warning?.Message, Does.Contain("N-terminal").And.Contain("UNIMOD:739"));
    }

    [Test]
    public void TryCleanSequence_NullSourceParser_UsesModelsOwnConverterParser()
    {
        var model = new KoinaModelHarness(KoinaModelHarness.BuildConverter(new HashSet<int> { 35 }), allowedUnimodIds: new HashSet<int> { 35 });

        var result = model.TryCleanWithParser("PEPM[Common Variable:Oxidation on M]IDE", null, out var apiSequence, out _);

        Assert.That(result, Is.Not.Null);
        Assert.That(apiSequence, Does.Contain("UNIMOD:35"));
    }

    [Test]
    public void TryCleanSequence_ExplicitSourceParser_OverridesConvertersOwnParser()
    {
        // The caller's parser supplies both the brackets that separate modifications from residues and the parse:
        // under the converter's own mzLib schema, "{note}" would be residues and fail.
        var fakeParser = new FakeSequenceParser(_ => CanonicalSequence.Unmodified("PEPTIDEK", "fake"), new BraceSchema());
        var model = new KoinaModelHarness(KoinaModelHarness.BuildConverter(new HashSet<int>()));

        var result = model.TryCleanWithParser("PEPT{note}IDEK", fakeParser, out var apiSequence, out var warning);

        Assert.That(result, Is.EqualTo("PEPTIDEK"));
        Assert.That(apiSequence, Is.EqualTo("PEPTIDEK"));
        Assert.That(warning, Is.Null);
    }

    [Test]
    public void TryCleanSequence_PreIdentifiedUnimodIdOutsideAllowList_ReturnNullRejectsBeforeSerialization()
    {
        // Simulates a ProForma-style "UNIMOD:N" token that already carries a resolved UnimodId.
        // UnimodSequenceSerializer.ShouldResolveMod skips lookup for such mods, so without the
        // pre-serialization allow-list check this would bypass AllowedUnimodIds and reach Koina.
        var preIdentified = CanonicalModification.AtResidue(3, 'M', "UNIMOD:35", unimodId: 35);
        var fakeParser = new FakeSequenceParser(_ =>
            CanonicalSequence.Unmodified("PEPMIDE", "fake").WithModification(preIdentified));
        var model = new KoinaModelHarness(
            KoinaModelHarness.BuildConverter(new HashSet<int> { 4 }),
            allowedUnimodIds: new HashSet<int> { 4 }); // 35 not allowed

        var result = model.TryCleanWithParser("PEPM[UNIMOD:35]IDE", fakeParser, out var apiSequence, out var warning);

        Assert.That(result, Is.Null);
        Assert.That(apiSequence, Is.Null);
        Assert.That(warning, Is.Not.Null);
        Assert.That(warning!.Message, Does.Contain("UNIMOD:35"));
    }

    [Test]
    public void TryCleanSequence_PreIdentifiedUnimodIdWithinAllowList_ReachesSerializer()
    {
        var preIdentified = CanonicalModification.AtResidue(3, 'M', "UNIMOD:35", unimodId: 35);
        var fakeParser = new FakeSequenceParser(_ =>
            CanonicalSequence.Unmodified("PEPMIDE", "fake").WithModification(preIdentified));
        var model = new KoinaModelHarness(
            KoinaModelHarness.BuildConverter(new HashSet<int> { 35 }),
            allowedUnimodIds: new HashSet<int> { 35 });

        var result = model.TryCleanWithParser("PEPM[UNIMOD:35]IDE", fakeParser, out var apiSequence, out var warning);

        Assert.That(result, Is.Not.Null);
        Assert.That(apiSequence, Does.Contain("UNIMOD:35"));
    }

    // Oxidation (UNIMOD:35) is outside the allow-list. The mzLib parser leaves UnimodId unset, so the
    // serializer judges it; ProForma pre-identifies it, so the explicit allow-list check does. Every
    // mode must end the same way for both sources.
    [TestCase(SequenceConversionHandlingMode.ReturnNull, null)]
    [TestCase(SequenceConversionHandlingMode.RemoveIncompatibleElements, "PEPMIDE")]
    [TestCase(SequenceConversionHandlingMode.UsePrimarySequence, "PEPMIDE")]
    public void TryCleanSequence_DisallowedModification_SameOutcomeFromMzLibAndProFormaSources(SequenceConversionHandlingMode mode, string? expected)
    {
        var model = new KoinaModelHarness(KoinaModelHarness.BuildConverter(new HashSet<int> { 4 }), mode, new HashSet<int> { 4 });

        var fromMzLib = model.TryCleanWithParser("PEPM[Common Variable:Oxidation on M]IDE", null, out var mzLibApi, out var mzLibWarning);
        var fromProForma = model.TryCleanWithParser("PEPM[UNIMOD:35]IDE", ProFormaSequenceParser.Instance, out var proFormaApi, out var proFormaWarning);

        Assert.That(fromMzLib, Is.EqualTo(expected));
        Assert.That(fromProForma, Is.EqualTo(expected));
        Assert.That(proFormaApi, Is.EqualTo(mzLibApi));
        Assert.That(mzLibWarning, Is.Not.Null);
        Assert.That(proFormaWarning, Is.Not.Null);
    }

    [Test]
    public void TryCleanSequence_DisallowedModificationThrowMode_ThrowsFromMzLibAndProFormaSources()
    {
        var model = new KoinaModelHarness(
            KoinaModelHarness.BuildConverter(new HashSet<int> { 4 }),
            SequenceConversionHandlingMode.ThrowException,
            new HashSet<int> { 4 });

        Assert.Throws<ArgumentException>(() => model.TryCleanWithParser("PEPM[Common Variable:Oxidation on M]IDE", null, out _, out _));
        Assert.Throws<ArgumentException>(() => model.TryCleanWithParser("PEPM[UNIMOD:35]IDE", ProFormaSequenceParser.Instance, out _, out _));
    }

    [Test]
    public void TryCleanSequence_RemoveIncompatibleElements_KeepsAllowedPreIdentifiedModification()
    {
        // Only the disallowed modification is dropped; an allowed one on the same peptide survives.
        var model = new KoinaModelHarness(
            KoinaModelHarness.BuildConverter(new HashSet<int> { 4 }),
            SequenceConversionHandlingMode.RemoveIncompatibleElements,
            new HashSet<int> { 4 });

        var result = model.TryCleanWithParser("PEPM[UNIMOD:35]IDEC[UNIMOD:4]", ProFormaSequenceParser.Instance, out _, out var warning);

        Assert.That(result, Is.EqualTo("PEPMIDEC[UNIMOD:4]"));
        Assert.That(warning?.Message, Does.Contain("UNIMOD:35"));
    }

    [Test]
    public void TryCleanSequence_InvalidBaseSequenceThrowMode_Throws()
    {
        var model = new KoinaModelHarness(
            KoinaModelHarness.BuildConverter(new HashSet<int>()),
            SequenceConversionHandlingMode.ThrowException);

        Assert.Throws<ArgumentException>(() => model.TryClean("PEP*TIDE", out _, out _));
    }

    [Test]
    public void TryCleanSequence_UsePrimarySequence_StripsModificationsAndWarns()
    {
        var model = new KoinaModelHarness(
            KoinaModelHarness.BuildConverter(new HashSet<int> { 35 }),
            SequenceConversionHandlingMode.UsePrimarySequence,
            new HashSet<int> { 35 });

        var result = model.TryClean("PEPM[Common Variable:Oxidation on M]IDE", out var apiSequence, out var warning);

        Assert.That(result, Is.Not.Null);
        Assert.That(apiSequence, Is.EqualTo(result));
        Assert.That(result, Does.Not.Contain("["));
        Assert.That(result, Is.EqualTo("PEPMIDE"));
        Assert.That(warning, Is.Not.Null);
        Assert.That(warning!.Message, Does.Contain("removed"));
    }

    [Test]
    public void TryCleanSequence_UnsupportedModification_ReturnsWarning()
    {
        var model = new KoinaModelHarness(KoinaModelHarness.BuildConverter(new HashSet<int>()));

        var result = model.TryClean("PEPM[Common Variable:Oxidation on M]IDE", out var apiSequence, out var warning);

        Assert.That(result, Is.Null);
        Assert.That(apiSequence, Is.Null);
        Assert.That(warning, Is.Not.Null);
        Assert.That(warning!.Message, Does.Contain("unsupported modifications"));
    }

    [Test]
    public void TryCleanSequence_WhenParseReturnsNull_BuildsWarningAndReturnsNull()
    {
        var converter = new FakeSequenceConverter(
            parse: _ => null,
            serialize: _ => "PEPTIDE");
        var model = new KoinaModelHarness(converter);

        var result = model.TryClean("PEPTIDE", out var apiSequence, out var warning);

        Assert.That(result, Is.Null);
        Assert.That(apiSequence, Is.Null);
        Assert.That(warning, Is.Null);
    }

    [Test]
    public void TryCleanSequence_WhenParseThrows_BuildsWarningAndReturnsNull()
    {
        var converter = new FakeSequenceConverter(
            parse: _ => throw new SequenceConversionException("parse failed", ConversionFailureReason.InvalidSequence),
            serialize: _ => "PEPTIDE");
        var model = new KoinaModelHarness(converter);

        var result = model.TryClean("PEPTIDE", out var apiSequence, out var warning);

        Assert.That(result, Is.Null);
        Assert.That(apiSequence, Is.Null);
        Assert.That(warning, Is.Not.Null);
        Assert.That(warning!.Message, Does.Contain("parse failed"));
    }

    [Test]
    public void TryCleanSequence_WhenSerializeThrows_BuildsWarningAndReturnsNull()
    {
        var converter = new FakeSequenceConverter(
            parse: _ => CanonicalSequence.Unmodified("PEPTIDE", "fake"),
            serialize: _ => throw new SequenceConversionException("serialize failed", ConversionFailureReason.InvalidSequence));
        var model = new KoinaModelHarness(converter);

        var result = model.TryClean("PEPTIDE", out var apiSequence, out var warning);

        Assert.That(result, Is.Null);
        Assert.That(apiSequence, Is.Null);
        Assert.That(warning, Is.Not.Null);
        Assert.That(warning!.Message, Does.Contain("serialize failed"));
    }

    [TestCase("PEPTIDE-[Unimod:Amidated on X]", "PEPTIDE-[UNIMOD:2]")]
    [TestCase("[Unimod:Acetyl on X]PEPTIDE-[Unimod:Amidated on X]", "[UNIMOD:1]PEPTIDE-[UNIMOD:2]")]
    public void TryCleanSequence_MzLibCTerminalModification_SurvivesTheRawBaseSequenceCheck(
        string mzLibSequence, string expectedApiSequence)
    {
        // The raw check used to strip only the brackets, leaving "PEPTIDE-", which the amino-acid pattern
        // rejects before the converter is ever asked. Under ReturnNull that failure was silent -- no
        // warning, so a caller could not tell it apart from a sequence that was never valid. Run against
        // the REAL converter, so this asserts the whole of TryCleanSequence and not just the pre-check.
        var model = new KoinaModelHarness(KoinaModelHarness.BuildAcceptAllConverter(), acceptsAllUnimodModifications: true);

        var result = model.TryClean(mzLibSequence, out var apiSequence, out var warning);

        Assert.That(result, Is.EqualTo(expectedApiSequence));
        Assert.That(apiSequence, Is.EqualTo(expectedApiSequence));
        Assert.That(warning, Is.Null);
    }

    [Test]
    public void TryCleanSequence_CTerminalModificationTheModelDoesNotAllow_FailsWithANamedWarning()
    {
        // Clearing the raw check is not the same as being predicted. A C-terminal group the model's
        // lookup does not know still fails -- but now it fails at the converter, which says WHICH
        // modification it was, instead of being dropped by the pre-check with no warning at all.
        var model = new KoinaModelHarness(KoinaModelHarness.BuildConverter(new HashSet<int> { 35 }));

        var result = model.TryClean("PEPTIDE-[Unimod:Amidated on X]", out var apiSequence, out var warning);

        Assert.That(result, Is.Null);
        Assert.That(apiSequence, Is.Null);
        Assert.That(warning, Is.Not.Null);
        Assert.That(warning!.Message, Does.Contain("Amidated"));
    }

    // A bare separator, a separator anywhere but the C-terminus, and ProForma shapes the Koina
    // converters cannot read. MzLibSequenceParser would parse
    // "[Unimod:Acetyl on X]-[Unimod:Amidated on X]PEPTIDE" as N-terminal acetylation PLUS a
    // C-terminal amidation -- a different peptide than was asked for, with no warning.
    // Being stopped here, before parsing, is the correct outcome for all of them.
    [TestCase("-PEPTIDE")]
    [TestCase("PEPTIDE-")]
    [TestCase("PEP-[UNIMOD:1]TIDE")]
    [TestCase("[UNIMOD:1]-PEP*TIDE")]
    [TestCase("[UNIMOD:1]-PEPTIDE")]
    [TestCase("[UNIMOD:1][UNIMOD:34]-PEPTIDE")]
    [TestCase("[Unimod:Acetyl on X]-[Unimod:Amidated on X]PEPTIDE")]
    public void TryCleanSequence_SeparatorThatIsNotAnMzLibCTerminalModification_RejectedBeforeParsing(
        string sequence)
    {
        bool parseCalled = false;
        var converter = new FakeSequenceConverter(
            parse: _ => { parseCalled = true; return CanonicalSequence.Unmodified("PEPTIDE", "fake"); },
            serialize: _ => "PEPTIDE");
        var model = new KoinaModelHarness(converter);

        var result = model.TryClean(sequence, out var apiSequence, out _);

        Assert.That(parseCalled, Is.False);
        Assert.That(result, Is.Null);
        Assert.That(apiSequence, Is.Null);
    }

    [Test]
    public void TryCleanSequence_AcceptAllConverter_SerializesKnownModification()
    {
        // CreateUnimodConverterAcceptAll backs ms2pip / AlphaPeptDeep, which accept any UNIMOD mod
        // regardless of the model's AllowedUnimodIds set; such models declare AcceptsAllUnimodModifications.
        var model = new KoinaModelHarness(KoinaModelHarness.BuildAcceptAllConverter(), acceptsAllUnimodModifications: true);

        var result = model.TryClean("PEPM[Common Variable:Oxidation on M]IDE", out var apiSequence, out _);

        Assert.That(result, Is.Not.Null);
        Assert.That(apiSequence, Does.Contain("UNIMOD:35"));
    }

    [Test]
    public void TryCleanSequence_ExceedsMaxLength_ReturnsNull()
    {
        var model = new KoinaModelHarness(KoinaModelHarness.BuildConverter(new HashSet<int>()));

        var result = model.TryClean(new string('A', 60), out var apiSequence, out _);

        Assert.That(result, Is.Null);
        Assert.That(apiSequence, Is.Null);
    }

    [Test]
    public void TryGetUnimodId_HandlesAccessionAndDatabaseReferenceBranches()
    {
        var method = typeof(KoinaModelBase<string, string>).GetMethod("TryGetUnimodId", BindingFlags.NonPublic | BindingFlags.Static)!;

        var byAccession = new Modification(_originalId: "x", _accession: "UNIMOD:35", _target: null);
        var byDbReference = new Modification(
            _originalId: "x",
            _accession: null,
            _target: null,
            _databaseReference: new Dictionary<string, IList<string>>
            {
                { "OTHER", new List<string> { "x" } },
                { "UNIMOD", new List<string> { "UNIMOD::4" } }
            });
        var noId = new Modification(
            _originalId: "x",
            _accession: null,
            _target: null,
            _databaseReference: new Dictionary<string, IList<string>>
            {
                { "UNIMOD", new List<string>() }
            });

        var args1 = new object[] { byAccession, 0 };
        var args2 = new object[] { byDbReference, 0 };
        var args3 = new object[] { noId, 0 };

        Assert.That((bool)method.Invoke(null, args1)!, Is.True);
        Assert.That((int)args1[1], Is.EqualTo(35));
        Assert.That((bool)method.Invoke(null, args2)!, Is.True);
        Assert.That((int)args2[1], Is.EqualTo(4));
        Assert.That((bool)method.Invoke(null, args3)!, Is.False);
        Assert.That((int)args3[1], Is.EqualTo(-1));
    }
}
