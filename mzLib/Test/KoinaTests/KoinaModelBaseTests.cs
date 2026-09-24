using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Reflection;
using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;

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
        }

        public string FormatName => "fake-fake";
        public string SourceFormatName => "fake";
        public string TargetFormatName => "fake";
        public ISequenceParser Parser => null!;
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

    private sealed class FakeSequenceParser : ISequenceParser
    {
        private readonly Func<string, CanonicalSequence?> _parse;

        public FakeSequenceParser(Func<string, CanonicalSequence?> parse)
        {
            _parse = parse;
        }

        public string FormatName => "fake";
        public SequenceFormatSchema Schema => MzLibSequenceFormatSchema.Instance;
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
            bool acceptsAllUnimodModifications = false)
            : base(converter)
        {
            ModHandlingMode = modHandlingMode;
            AllowedUnimodIds = allowedUnimodIds ?? new HashSet<int>();
            AcceptsAllUnimodModifications = acceptsAllUnimodModifications;
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
        // MzLibSequenceParser drops '*' with a warning rather than failing the parse, so
        // TryCleanSequence promotes that warning to an actionable rejection instead of silently
        // treating "PEP*TIDE" as "PEPTIDE".
        Assert.That(warning, Is.Not.Null);
        Assert.That(warning!.Message, Does.Contain("unsupported or ignored syntax"));
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
        // A fake parser that returns a fixed CanonicalSequence proves the caller-supplied
        // SequenceParser is actually consulted instead of the converter's own MzLib parser.
        var fakeParser = new FakeSequenceParser(_ => CanonicalSequence.Unmodified("PEPTIDEK", "fake"));
        var model = new KoinaModelHarness(KoinaModelHarness.BuildConverter(new HashSet<int>()));

        var result = model.TryCleanWithParser("this mzLib parser would reject this string", fakeParser, out var apiSequence, out var warning);

        Assert.That(result, Is.EqualTo("PEPTIDEK"));
        Assert.That(apiSequence, Is.EqualTo("PEPTIDEK"));
        Assert.That(warning, Is.Null);
    }

    [Test]
    public void TryCleanSequence_PreIdentifiedUnimodIdOutsideAllowList_IsRejectedBeforeSerialization()
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

    [Test]
    public void TryCleanSequence_AcceptsAllUnimodModifications_SkipsAllowListCheckForPreIdentifiedId()
    {
        var preIdentified = CanonicalModification.AtResidue(3, 'M', "UNIMOD:9999", unimodId: 9999); // not a real id
        var fakeParser = new FakeSequenceParser(_ =>
            CanonicalSequence.Unmodified("PEPMIDE", "fake").WithModification(preIdentified));
        var model = new KoinaModelHarness(
            KoinaModelHarness.BuildAcceptAllConverter(),
            acceptsAllUnimodModifications: true);

        // AcceptsAllUnimodModifications skips the pre-serialization allow-list check, so no early
        // rejection naming UNIMOD:9999 is produced here (serialization-time lookup is not this test's concern).
        var result = model.TryCleanWithParser("PEPM[UNIMOD:9999]IDE", fakeParser, out _, out var warning);

        Assert.That(warning is null || !warning.Message.Contains("Sequence contains unsupported modification(s): UNIMOD:9999"));
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

    [Test]
    public void TryCleanSequence_AcceptAllConverter_SerializesKnownModification()
    {
        // CreateUnimodConverterAcceptAll backs ms2pip / AlphaPeptDeep, which accept any UNIMOD mod
        // regardless of the model's AllowedUnimodIds set.
        var model = new KoinaModelHarness(KoinaModelHarness.BuildAcceptAllConverter());

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
