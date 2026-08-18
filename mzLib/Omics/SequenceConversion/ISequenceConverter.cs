namespace Omics.SequenceConversion;

/// <summary>
/// Wraps a parser/serializer pair for a source-to-target conversion.
/// Provides a unified API for parse/serialize/convert operations.
/// </summary>
public interface ISequenceConverter
{
    /// <summary>
    /// Gets the unique name of the converter in the format "{sourceFormat}-{targetFormat}".
    /// </summary>
    string FormatName { get; }

    /// <summary>
    /// Gets the source format handled by this converter.
    /// </summary>
    string SourceFormatName { get; }

    /// <summary>
    /// Gets the target format handled by this converter.
    /// </summary>
    string TargetFormatName { get; }

    /// <summary>
    /// Gets the parser for the source format.
    /// </summary>
    ISequenceParser Parser { get; }

    /// <summary>
    /// Gets the serializer for the target format.
    /// </summary>
    ISequenceSerializer Serializer { get; }


    /// <summary>
    /// Parses an input string into a canonical sequence.
    /// </summary>
    CanonicalSequence? Parse(
        string input,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException);

    /// <summary>
    /// Serializes a canonical sequence into this format.
    /// </summary>
    string? Serialize(
        CanonicalSequence sequence,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException);

    /// <summary>
    /// Converts an input string in this format to a canonical sequence and back.
    /// </summary>
    string? Convert(
        string input,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException);
}
