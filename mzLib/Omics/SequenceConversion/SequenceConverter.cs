namespace Omics.SequenceConversion;

/// <summary>
/// Default implementation of <see cref="ISequenceConverter"/>.
/// </summary>
public sealed class SequenceConverter : ISequenceConverter
{
    public SequenceConverter(ISequenceParser? parser, ISequenceSerializer? serializer)
    {
        if (parser == null && serializer == null)
        {
            throw new ArgumentException("A sequence converter requires a parser or serializer.", nameof(parser));
        }

        Parser = parser;
        Serializer = serializer;

        var formatName = serializer?.FormatName ?? parser!.FormatName;
        if (parser != null && serializer != null &&
            !string.Equals(parser.FormatName, serializer.FormatName, StringComparison.OrdinalIgnoreCase))
        {
            throw new ArgumentException(
                $"Parser format '{parser.FormatName}' does not match serializer format '{serializer.FormatName}'.",
                nameof(serializer));
        }

        FormatName = formatName;
    }

    public string FormatName { get; }

    public ISequenceParser? Parser { get; }

    public ISequenceSerializer? Serializer { get; }

    public bool CanParse => Parser != null;

    public bool CanSerialize => Serializer != null;

    public CanonicalSequence? Parse(
        string input,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
    {
        warnings ??= new ConversionWarnings();

        if (Parser == null)
        {
            return SequenceConversionHelpers.HandleParserError(
                warnings,
                mode,
                ConversionFailureReason.UnknownFormat,
                $"No parser registered for format '{FormatName}'.");
        }

        return Parser.Parse(input, warnings, mode);
    }

    public string? Serialize(
        CanonicalSequence sequence,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
    {
        warnings ??= new ConversionWarnings();

        if (Serializer == null)
        {
            return SequenceConversionHelpers.HandleSerializerError(
                warnings,
                mode,
                ConversionFailureReason.UnknownFormat,
                $"No serializer registered for format '{FormatName}'.");
        }

        return Serializer.Serialize(sequence, warnings, mode);
    }

    public string? Convert(
        string input,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
    {
        warnings ??= new ConversionWarnings();

        var canonical = Parse(input, warnings, mode);
        if (canonical == null)
        {
            return null;
        }

        return Serialize(canonical.Value, warnings, mode);
    }
}
