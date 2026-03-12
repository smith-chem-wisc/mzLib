namespace Omics.SequenceConversion;

/// <summary>
/// Wraps a parser/serializer pair for a single format.
/// Provides a unified API for parse/serialize/convert operations.
/// </summary>
public interface ISequenceConverter
{
    /// <summary>
    /// Gets the unique name of the format handled by this converter.
    /// </summary>
    string FormatName { get; }

    /// <summary>
    /// Gets the parser for this format (null if parsing is not supported).
    /// </summary>
    ISequenceParser? Parser { get; }

    /// <summary>
    /// Gets the serializer for this format (null if serialization is not supported).
    /// </summary>
    ISequenceSerializer? Serializer { get; }

    /// <summary>
    /// Returns true if this converter can parse input strings.
    /// </summary>
    bool CanParse { get; }

    /// <summary>
    /// Returns true if this converter can serialize canonical sequences.
    /// </summary>
    bool CanSerialize { get; }

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
