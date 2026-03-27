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

    private sealed class KoinaModelHarness : KoinaModelBase<string, string>
    {
        public KoinaModelHarness(
            ISequenceConverter converter,
            SequenceConversionHandlingMode modHandlingMode = SequenceConversionHandlingMode.ReturnNull,
            IReadOnlySet<int>? allowedUnimodIds = null)
            : base(converter)
        {
            ModHandlingMode = modHandlingMode;
            AllowedUnimodIds = allowedUnimodIds ?? new HashSet<int>();
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

        protected override List<Dictionary<string, object>> ToBatchedRequests(List<string> validInputs)
        {
            return new List<Dictionary<string, object>>();
        }

        public string? TryClean(string sequence, out string? apiSequence, out WarningException? warning)
        {
            return TryCleanSequence(sequence, out apiSequence, out warning);
        }

        public static ISequenceConverter BuildConverter(IReadOnlySet<int> allowedUnimodIds)
        {
            return CreateUnimodConverter(UnimodSequenceFormatSchema.Instance, allowedUnimodIds);
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
        Assert.That(warning, Is.Null);
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
