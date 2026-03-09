namespace Omics.SequenceConversion;

/// <summary>
/// Base class for defining the syntax rules for a sequence format.
/// Used by parsers and serializers to understand how sequences are structured.
/// Derive from this class to create format-specific schemas with custom behavior.
/// </summary>
public abstract class SequenceFormatSchema
{
    /// <summary>
    /// Unique identifier for this format (e.g., "mzLib", "Chronologer", "UNIMOD").
    /// </summary>
    public abstract string FormatName { get; }

    /// <summary>
    /// Character that opens a modification annotation (e.g., '[' or '(').
    /// </summary>
    public abstract char ModOpenBracket { get; }

    /// <summary>
    /// Character that closes a modification annotation (e.g., ']' or ')').
    /// </summary>
    public abstract char ModCloseBracket { get; }

    /// <summary>
    /// Optional separator after N-terminal modification (e.g., "-").
    /// Null if format doesn't use N-terminal separators.
    /// Empty string if N-terminal mods directly precede the sequence without separator.
    /// </summary>
    public abstract string? NTermSeparator { get; }

    /// <summary>
    /// Optional separator before C-terminal modification (e.g., "-").
    /// Null if format doesn't use C-terminal separators.
    /// </summary>
    public abstract string? CTermSeparator { get; }

    /// <summary>
    /// Returns true if this format supports N-terminal modification annotations.
    /// </summary>
    public bool SupportsNTerminalMods => NTermSeparator != null;

    /// <summary>
    /// Returns true if this format supports C-terminal modification annotations.
    /// </summary>
    public bool SupportsCTerminalMods => CTermSeparator != null;

    /// <summary>
    /// Checks if a character is a modification bracket (open or close).
    /// </summary>
    public bool IsModBracket(char c) => c == ModOpenBracket || c == ModCloseBracket;

    public override string ToString() => FormatName;
}
