namespace Omics.SequenceConversion;

/// <summary>
/// Parses a string representation of a sequence into a <see cref="CanonicalSequence"/>.
/// Each implementation handles a specific input format (e.g., mzLib, UNIMOD, mass-shift).
/// </summary>
public interface ISequenceParser
{
    /// <summary>
    /// Gets the unique name of the format this parser handles (e.g., "mzLib", "UNIMOD").
    /// </summary>
    string FormatName { get; }

    /// <summary>
    /// Gets the schema describing the syntax of this format.
    /// </summary>
    SequenceFormatSchema Schema { get; }

    /// <summary>
    /// Determines whether this parser can likely parse the given input.
    /// This is a fast heuristic check, not a full validation.
    /// </summary>
    /// <param name="input">The input string to check.</param>
    /// <returns>True if this parser can likely handle the input; otherwise, false.</returns>
    bool CanParse(string input);

    /// <summary>
    /// Parses the input string into a <see cref="CanonicalSequence"/>.
    /// </summary>
    /// <param name="input">The input string to parse.</param>
    /// <param name="warnings">Optional. Accumulates warnings and errors during parsing.
    /// If null, a new instance will be created internally if issues occur.</param>
    /// <param name="mode">Specifies how to handle parsing issues (throw, return null, etc.).</param>
    /// <returns>The parsed sequence, or null if parsing failed and mode allows null returns.</returns>
    /// <exception cref="SequenceConversionException">Thrown when parsing fails and mode is ThrowException.</exception>
    CanonicalSequence? Parse(
        string input,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException);
}
