namespace Omics.SequenceConversion;

/// <summary>
/// Orchestrates sequence conversions between different formats.
/// Manages registered parsers and serializers, and provides a unified API for conversions.
/// </summary>
public interface ISequenceConversionService
{
    /// <summary>
    /// Gets the names of all registered source formats (formats that can be parsed).
    /// </summary>
    IReadOnlyCollection<string> AvailableSourceFormats { get; }

    /// <summary>
    /// Gets the names of all registered target formats (formats that can be serialized to).
    /// </summary>
    IReadOnlyCollection<string> AvailableTargetFormats { get; }

    /// <summary>
    /// Gets the names of all registered converters (formats with parser/serializer wrappers).
    /// </summary>
    IReadOnlyCollection<string> AvailableConverters { get; }

    /// <summary>
    /// Parses an input string from the specified source format into a <see cref="CanonicalSequence"/>.
    /// </summary>
    /// <param name="input">The input string to parse.</param>
    /// <param name="sourceFormat">The name of the source format (e.g., "mzLib", "UNIMOD").</param>
    /// <param name="warnings">Optional. Accumulates warnings and errors during parsing.</param>
    /// <param name="mode">Specifies how to handle parsing issues.</param>
    /// <returns>The parsed sequence, or null if parsing failed and mode allows null returns.</returns>
    CanonicalSequence? Parse(
        string input,
        string sourceFormat,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException);

    /// <summary>
    /// Serializes a <see cref="CanonicalSequence"/> to the specified target format.
    /// </summary>
    /// <param name="sequence">The sequence to serialize.</param>
    /// <param name="targetFormat">The name of the target format (e.g., "Chronologer", "MassShift").</param>
    /// <param name="warnings">Optional. Accumulates warnings and errors during serialization.</param>
    /// <param name="mode">Specifies how to handle serialization issues.</param>
    /// <returns>The serialized string, or null if serialization failed and mode allows null returns.</returns>
    string? Serialize(
        CanonicalSequence sequence,
        string targetFormat,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException);

    /// <summary>
    /// Converts an input string from one format to another in a single operation.
    /// This is equivalent to calling Parse followed by Serialize.
    /// </summary>
    /// <param name="input">The input string to convert.</param>
    /// <param name="sourceFormat">The name of the source format.</param>
    /// <param name="targetFormat">The name of the target format.</param>
    /// <param name="warnings">Optional. Accumulates warnings and errors during conversion.</param>
    /// <param name="mode">Specifies how to handle conversion issues.</param>
    /// <returns>The converted string, or null if conversion failed and mode allows null returns.</returns>
    string? Convert(
        string input,
        string sourceFormat,
        string targetFormat,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException);

    /// <summary>
    /// Attempts to auto-detect the format of the input string.
    /// </summary>
    /// <param name="input">The input string to analyze.</param>
    /// <returns>The detected format name, or null if format cannot be determined.</returns>
    string? DetectFormat(string input);

    /// <summary>
    /// Registers a parser for a specific format.
    /// </summary>
    /// <param name="parser">The parser to register.</param>
    void RegisterParser(ISequenceParser parser);

    /// <summary>
    /// Registers a serializer for a specific format.
    /// </summary>
    /// <param name="serializer">The serializer to register.</param>
    void RegisterSerializer(ISequenceSerializer serializer);

    /// <summary>
    /// Registers a converter for a specific format.
    /// </summary>
    /// <param name="converter">The converter to register.</param>
    void RegisterConverter(ISequenceConverter converter);

    /// <summary>
    /// Gets the converter for a specific format, or null if not registered.
    /// </summary>
    ISequenceConverter? GetConverter(string formatName);
}
