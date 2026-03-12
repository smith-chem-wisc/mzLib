namespace Omics.SequenceConversion;

/// <summary>
/// Default implementation of <see cref="ISequenceConverter"/>.
/// </summary>
public sealed class SequenceConverter : ISequenceConverter
{
    public SequenceConverter(ISequenceParser parser, ISequenceSerializer serializer)
    {
        ArgumentNullException.ThrowIfNull(parser);
        ArgumentNullException.ThrowIfNull(serializer);

        Parser = parser;
        Serializer = serializer;

        SourceFormatName = parser.FormatName;
        TargetFormatName = serializer.FormatName;
        FormatName = $"{SourceFormatName}-{TargetFormatName}";
    }

    public string FormatName { get; }

    public string SourceFormatName { get; }

    public string TargetFormatName { get; }

    public ISequenceParser Parser { get; }

    public ISequenceSerializer Serializer { get; }


    public CanonicalSequence? Parse(
        string input,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
    {
        warnings ??= new ConversionWarnings();

        return Parser.Parse(input, warnings, mode);
    }

    public string? Serialize(
        CanonicalSequence sequence,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
    {
        warnings ??= new ConversionWarnings();

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
