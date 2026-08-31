namespace Omics.SequenceConversion;

/// <summary>
/// Base class for defining the syntax rules for a sequence format.
/// Used by parsers and serializers to understand how sequences are structured.
/// Derive from this class to create format-specific schemas with custom behavior.
/// </summary>
public abstract class SequenceFormatSchema(char modOpen = '[', char modClosed = ']', string? nTermSeparator = "", string? cTermSeparator = "-")
{
    /// <summary>
    /// Unique identifier for this format (e.g., "mzLib", "Chronologer", "UNIMOD").
    /// </summary>
    public abstract string FormatName { get; }

    /// <summary>
    /// Character that opens a modification annotation (e.g., '[' or '(').
    /// </summary>
    public virtual char ModOpenBracket { get; init; } = modOpen;

    /// <summary>
    /// Character that closes a modification annotation (e.g., ']' or ')').
    /// </summary>
    public virtual char ModCloseBracket { get; init; } = modClosed;

    /// <summary>
    /// Optional separator after N-terminal modification (e.g., "-").
    /// Null if format doesn't use N-terminal separators.
    /// Empty string if N-terminal mods directly precede the sequence without separator. 
    /// <remarks>
    /// Empty string indicates N-terminal modifications directly precede the sequence with no separator.
    /// Example: "[Acetyl]PEPTIDE" not "[Acetyl]-PEPTIDE"
    /// </remarks>
    /// </summary>
    public virtual string? NTermSeparator { get; init; } = nTermSeparator;

    /// <summary>
    /// Optional separator before C-terminal modification (e.g., "-").
    /// Null if format doesn't use C-terminal separators.
    /// </summary>
    public virtual string? CTermSeparator { get; init; } = cTermSeparator;

    /// <summary>
    /// Returns true if this format supports N-terminal modification annotations.
    /// </summary>
    public bool SupportsNTerminalMods => NTermSeparator != null;

    /// <summary>
    /// Returns true if this format supports C-terminal modification annotations.
    /// </summary>
    public bool SupportsCTerminalMods => CTermSeparator != null;

    public override string ToString() => FormatName;
}
