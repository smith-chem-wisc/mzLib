namespace Omics.SequenceConversion;

/// <summary>
/// Implementation of <see cref="ISequenceConversionService"/> that manages sequence format conversions.
/// 
/// This service provides a centralized way to:
/// - Parse sequences from various input formats into <see cref="CanonicalSequence"/>
/// - Serialize <see cref="CanonicalSequence"/> to various output formats
/// - Convert directly between formats in a single operation
/// - Auto-detect input formats
/// 
/// The default instance comes pre-registered with mzLib and Chronologer formats.
/// Additional formats can be registered at runtime.
/// </summary>
public class SequenceConversionService : ISequenceConversionService
{
    private readonly Dictionary<string, ISequenceParser> _parsers = new(StringComparer.OrdinalIgnoreCase);
    private readonly Dictionary<string, ISequenceSerializer> _serializers = new(StringComparer.OrdinalIgnoreCase);
    private readonly Dictionary<string, ISequenceConverter> _converters = new(StringComparer.OrdinalIgnoreCase);

    /// <summary>
    /// Default instance with mzLib and Chronologer formats pre-registered.
    /// </summary>
    public static SequenceConversionService Default { get; } = CreateDefault();

    /// <summary>
    /// Creates a new empty SequenceConversionService.
    /// Use <see cref="Default"/> for a pre-configured instance.
    /// </summary>
    public SequenceConversionService() { }

    /// <summary>
    /// Creates the default service with standard parsers and serializers registered.
    /// </summary>
    private static SequenceConversionService CreateDefault()
    {
        var service = new SequenceConversionService();

        // Register mzLib format
        service.RegisterParser(MzLibSequenceParser.Instance);
        service.RegisterSerializer(MzLibSequenceSerializer.Instance);
        service.RegisterConverter(new SequenceConverter(MzLibSequenceParser.Instance, MzLibSequenceSerializer.Instance));

        // Register mass shift format
        service.RegisterParser(MassShiftSequenceParser.Instance);
        service.RegisterSerializer(MassShiftSequenceSerializer.Instance);
        service.RegisterConverter(new SequenceConverter(MzLibSequenceParser.Instance, MassShiftSequenceSerializer.Instance));

        // Register Chronologer format (serializer only)
        service.RegisterSerializer(ChronologerSequenceSerializer.Instance);
        service.RegisterConverter(new SequenceConverter(MzLibSequenceParser.Instance, ChronologerSequenceSerializer.Instance));

        // Register Unimod format (serializer only)
        service.RegisterSerializer(UnimodSequenceSerializer.Instance);
        service.RegisterConverter(new SequenceConverter(MzLibSequenceParser.Instance, UnimodSequenceSerializer.Instance));

        // Register UniProt format (serializer only)
        service.RegisterSerializer(UniProtSequenceSerializer.Instance);
        service.RegisterConverter(new SequenceConverter(MzLibSequenceParser.Instance, UniProtSequenceSerializer.Instance));

        // Register EssentialSequence (serializer only with default (from MM) mod allowances are w)
        service.RegisterSerializer(EssentialSequenceSerializer.Instance);
        service.RegisterConverter(new SequenceConverter(MzLibSequenceParser.Instance, EssentialSequenceSerializer.Instance));

        return service;
    }

    /// <inheritdoc />
    public IReadOnlyCollection<string> AvailableSourceFormats => _parsers.Keys.ToList().AsReadOnly();

    /// <inheritdoc />
    public IReadOnlyCollection<string> AvailableTargetFormats => _serializers.Keys.ToList().AsReadOnly();

    /// <inheritdoc />
    public IReadOnlyCollection<string> AvailableConverters => BuildAvailableConverters();

    /// <inheritdoc />
    public void RegisterParser(ISequenceParser parser)
    {
        ArgumentNullException.ThrowIfNull(parser);
        _parsers[parser.FormatName] = parser;
    }

    /// <inheritdoc />
    public void RegisterSerializer(ISequenceSerializer serializer)
    {
        ArgumentNullException.ThrowIfNull(serializer);
        _serializers[serializer.FormatName] = serializer;
    }

    /// <inheritdoc />
    public void RegisterConverter(ISequenceConverter converter)
    {
        ArgumentNullException.ThrowIfNull(converter);
        _converters[converter.FormatName] = converter;
    }

    /// <inheritdoc />
    public CanonicalSequence? Parse(
        string input,
        string sourceFormat,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
    {
        warnings ??= new ConversionWarnings();

        if (string.IsNullOrWhiteSpace(input))
        {
            return HandleParseError(warnings, mode, ConversionFailureReason.InvalidSequence,
                "Input sequence is null or empty.");
        }

        if (!_parsers.TryGetValue(sourceFormat, out var parser))
        {
            return HandleParseError(warnings, mode, ConversionFailureReason.UnknownFormat,
                $"No parser registered for format '{sourceFormat}'. Available formats: {string.Join(", ", _parsers.Keys)}");
        }

        return parser.Parse(input, warnings, mode);
    }

    /// <inheritdoc />
    public string? Serialize(
        CanonicalSequence sequence,
        string targetFormat,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
    {
        warnings ??= new ConversionWarnings();

        if (string.IsNullOrEmpty(sequence.BaseSequence))
        {
            return HandleSerializeError(warnings, mode, ConversionFailureReason.InvalidSequence,
                "Sequence has no base sequence.");
        }

        if (!_serializers.TryGetValue(targetFormat, out var serializer))
        {
            return HandleSerializeError(warnings, mode, ConversionFailureReason.UnknownFormat,
                $"No serializer registered for format '{targetFormat}'. Available formats: {string.Join(", ", _serializers.Keys)}");
        }

        return serializer.Serialize(sequence, warnings, mode);
    }

    /// <inheritdoc />
    public string? Convert(
        string input,
        string sourceFormat,
        string targetFormat,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
    {
        warnings ??= new ConversionWarnings();

        var converter = GetConverter(sourceFormat, targetFormat);
        if (converter != null)
        {
            return converter.Convert(input, warnings, mode);
        }

        // Parse to canonical form
        var canonical = Parse(input, sourceFormat, warnings, mode);
        if (canonical == null)
            return null;

        // Serialize to target format
        return Serialize(canonical.Value, targetFormat, warnings, mode);
    }

    /// <inheritdoc />
    public string? DetectFormat(string input)
    {
        if (string.IsNullOrWhiteSpace(input))
            return null;

        // Try each parser's CanParse method to detect format
        // Return the first one that matches
        foreach (var parser in _parsers.Values)
        {
            if (parser.CanParse(input))
                return parser.FormatName;
        }

        return null;
    }

    /// <summary>
    /// Parses a sequence with auto-detection of the source format.
    /// </summary>
    /// <param name="input">The input string to parse.</param>
    /// <param name="warnings">Optional. Accumulates warnings and errors during parsing.</param>
    /// <param name="mode">Specifies how to handle parsing issues.</param>
    /// <returns>The parsed sequence, or null if parsing failed.</returns>
    public CanonicalSequence? ParseAutoDetect(
        string input,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
    {
        warnings ??= new ConversionWarnings();

        var detectedFormat = DetectFormat(input);
        if (detectedFormat == null)
        {
            return HandleParseError(warnings, mode, ConversionFailureReason.UnknownFormat,
                $"Could not auto-detect format for input: {(input.Length > 50 ? input.Substring(0, 50) + "..." : input)}");
        }

        return Parse(input, detectedFormat, warnings, mode);
    }

    /// <summary>
    /// Converts a sequence to the target format with auto-detection of the source format.
    /// </summary>
    /// <param name="input">The input string to convert.</param>
    /// <param name="targetFormat">The name of the target format.</param>
    /// <param name="warnings">Optional. Accumulates warnings and errors during conversion.</param>
    /// <param name="mode">Specifies how to handle conversion issues.</param>
    /// <returns>The converted string, or null if conversion failed.</returns>
    public string? ConvertAutoDetect(
        string input,
        string targetFormat,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
    {
        warnings ??= new ConversionWarnings();

        var canonical = ParseAutoDetect(input, warnings, mode);
        if (canonical == null)
            return null;

        return Serialize(canonical.Value, targetFormat, warnings, mode);
    }

    /// <summary>
    /// Gets the parser for a specific format, or null if not registered.
    /// </summary>
    public ISequenceParser? GetParser(string formatName)
    {
        return CanParseFormat(formatName) 
            ? _parsers.GetValueOrDefault(formatName) 
            : null;
    }

    /// <summary>
    /// Gets the serializer for a specific format, or null if not registered.
    /// </summary>
    public ISequenceSerializer? GetSerializer(string formatName)
    {
        return CanSerializeFormat(formatName) 
            ? _serializers.GetValueOrDefault(formatName) 
            : null;
    }

    /// <summary>
    /// Gets the converter for a specific source and target format, or null if unavailable.
    /// </summary>
    public ISequenceConverter? GetConverter(string sourceFormat, string targetFormat)
    {
        if (string.IsNullOrWhiteSpace(sourceFormat) || string.IsNullOrWhiteSpace(targetFormat))
        {
            return null;
        }

        var key = BuildConverterKey(sourceFormat, targetFormat);
        if (_converters.TryGetValue(key, out var converter))
        {
            return converter;
        }

        if (CanParseFormat(sourceFormat) && CanSerializeFormat(targetFormat))
        {
            converter = new SequenceConverter(_parsers[sourceFormat], _serializers[targetFormat]);
            _converters[key] = converter;
            return converter;
        }

        return null;
    }

    /// <summary>
    /// Checks if a specific source format is available for parsing.
    /// </summary>
    public bool CanParseFormat(string formatName) => _parsers.ContainsKey(formatName);

    /// <summary>
    /// Checks if a specific target format is available for serialization.
    /// </summary>
    public bool CanSerializeFormat(string formatName) => _serializers.ContainsKey(formatName);

    private static string BuildConverterKey(string sourceFormat, string targetFormat)
        => $"{sourceFormat}-{targetFormat}";

    private IReadOnlyCollection<string> BuildAvailableConverters()
    {
        if (_parsers.Count == 0 || _serializers.Count == 0)
        {
            return Array.Empty<string>();
        }

        var converters = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        foreach (var source in _parsers.Keys)
        {
            foreach (var target in _serializers.Keys)
            {
                converters.Add(BuildConverterKey(source, target));
            }
        }

        return converters.ToList().AsReadOnly();
    }

    #region Error Handling

    private static CanonicalSequence? HandleParseError(
        ConversionWarnings warnings,
        SequenceConversionHandlingMode mode,
        ConversionFailureReason reason,
        string message)
    {
        warnings.SetFailure(reason, message);

        return mode switch
        {
            SequenceConversionHandlingMode.ThrowException =>
                throw warnings.ToException(message),
            _ => null
        };
    }

    private static string? HandleSerializeError(
        ConversionWarnings warnings,
        SequenceConversionHandlingMode mode,
        ConversionFailureReason reason,
        string message)
    {
        warnings.SetFailure(reason, message);

        return mode switch
        {
            SequenceConversionHandlingMode.ThrowException =>
                throw warnings.ToException(message),
            _ => null
        };
    }

    #endregion
}
