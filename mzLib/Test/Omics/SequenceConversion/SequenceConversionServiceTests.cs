using NUnit.Framework;
using Omics.SequenceConversion;
using System;
using System.Collections.Generic;

namespace Test.Omics.SequenceConversion;

[TestFixture]
public class SequenceConversionServiceTests
{
    [Test]
    public void Parse_WithUnknownFormat_RecordsFailureAndReturnsNull()
    {
        var service = new SequenceConversionService();
        var warnings = new ConversionWarnings();

        var result = service.Parse("PEPTIDE", "does-not-exist", warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(result, Is.Null);
        Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.UnknownFormat));
        Assert.That(warnings.HasErrors, Is.True);
    }

    [Test]
    public void Convert_UsesRegisteredConverterWhenAvailable()
    {
        var parser = new StubParser("source", _ => true);
        var serializer = new StubSerializer("target");
        var converter = new StubConverter(parser, serializer, "converted");
        var service = new SequenceConversionService();
        service.RegisterConverter(converter);

        var result = service.Convert("input", parser.FormatName, serializer.FormatName, null, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(result, Is.EqualTo("converted:input"));
        Assert.That(converter.ConvertCallCount, Is.EqualTo(1));
    }

    [Test]
    public void GetConverter_BuildsHybridConverterAndCachesIt()
    {
        var parser = new StubParser("source", _ => true);
        var serializer = new StubSerializer("target");
        var service = new SequenceConversionService();
        service.RegisterParser(parser);
        service.RegisterSerializer(serializer);

        var converterOne = service.GetConverter("source", "target");
        var converterTwo = service.GetConverter("source", "target");

        Assert.That(converterOne, Is.Not.Null);
        Assert.That(converterTwo, Is.SameAs(converterOne));
        Assert.That(converterOne!.FormatName, Is.EqualTo("source-target"));
    }

    [Test]
    public void ParseAutoDetect_WhenNothingMatches_ReturnsNullAndFailure()
    {
        var service = new SequenceConversionService();
        var warnings = new ConversionWarnings();

        var result = service.ParseAutoDetect("???", warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(result, Is.Null);
        Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.UnknownFormat));
    }

    [Test]
    public void ConvertAutoDetect_UsesDetectedParserAndSerializer()
    {
        var parser = new StubParser("auto", input => input.StartsWith("AUTO:"), input =>
            CanonicalSequence.Unmodified(input.Split(':', 2)[1], "auto"));
        var serializer = new StubSerializer("target");
        var service = new SequenceConversionService();
        service.RegisterParser(parser);
        service.RegisterSerializer(serializer);

        var warnings = new ConversionWarnings();
        var result = service.ConvertAutoDetect("AUTO:SEQUENCE", serializer.FormatName, warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(result, Is.EqualTo("target:SEQUENCE"));
        Assert.That(warnings.HasFatalError, Is.False);
    }

    private sealed class StubParser : ISequenceParser
    {
        private readonly Func<string, bool> _canParse;
        private readonly Func<string, CanonicalSequence> _parse;

        public StubParser(string formatName, Func<string, bool> canParse)
            : this(formatName, canParse, input => CanonicalSequence.Unmodified(input, formatName))
        {
        }

        public StubParser(string formatName, Func<string, bool> canParse, Func<string, CanonicalSequence> parse)
        {
            FormatName = formatName;
            Schema = new StubSchema(formatName);
            _canParse = canParse;
            _parse = parse;
        }

        public string FormatName { get; }
        public SequenceFormatSchema Schema { get; }

        public bool CanParse(string input) => _canParse(input);

        public CanonicalSequence? Parse(
            string input,
            ConversionWarnings? warnings = null,
            SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException) => _parse(input);
    }

    private sealed class StubSerializer : ISequenceSerializer
    {
        public StubSerializer(string formatName, SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException)
        {
            FormatName = formatName;
            Schema = new StubSchema(formatName);
            HandlingMode = handlingMode;
        }

        public string FormatName { get; }
        public SequenceFormatSchema Schema { get; }
        public IModificationLookup? ModificationLookup => null;
        public SequenceConversionHandlingMode HandlingMode { get; }

        public bool CanSerialize(CanonicalSequence sequence) => true;
        public bool ShouldResolveMod(CanonicalModification mod) => false;

        public string? Serialize(
            CanonicalSequence sequence,
            ConversionWarnings? warnings = null,
            SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
            => $"{FormatName}:{sequence.BaseSequence}";
    }

    private sealed class StubConverter : ISequenceConverter
    {
        private readonly string _prefix;

        public StubConverter(ISequenceParser parser, ISequenceSerializer serializer, string prefix)
        {
            Parser = parser;
            Serializer = serializer;
            SourceFormatName = parser.FormatName;
            TargetFormatName = serializer.FormatName;
            FormatName = $"{SourceFormatName}-{TargetFormatName}";
            _prefix = prefix;
        }

        public string FormatName { get; }
        public string SourceFormatName { get; }
        public string TargetFormatName { get; }
        public ISequenceParser Parser { get; }
        public ISequenceSerializer Serializer { get; }
        public int ConvertCallCount { get; private set; }

        public CanonicalSequence? Parse(
            string input,
            ConversionWarnings? warnings = null,
            SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException) => Parser.Parse(input, warnings, mode);

        public string? Serialize(
            CanonicalSequence sequence,
            ConversionWarnings? warnings = null,
            SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException) => Serializer.Serialize(sequence, warnings, mode);

        public string? Convert(
            string input,
            ConversionWarnings? warnings = null,
            SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
        {
            ConvertCallCount++;
            return $"{_prefix}:{input}";
        }
    }

    private sealed class StubSchema : SequenceFormatSchema
    {
        public StubSchema(string formatName)
        {
            FormatName = formatName;
        }

        public override string FormatName { get; }
    }
}
