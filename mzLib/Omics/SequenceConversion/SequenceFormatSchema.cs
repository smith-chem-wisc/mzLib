namespace Omics.SequenceConversion;

/// <summary>
/// Defines the syntax rules for a sequence format.
/// Used by parsers and serializers to understand how sequences are structured.
/// </summary>
/// <param name="FormatName">Unique identifier for this format (e.g., "mzLib", "Chronologer", "UNIMOD").</param>
/// <param name="ModOpenBracket">Character that opens a modification annotation (e.g., '[' or '(').</param>
/// <param name="ModCloseBracket">Character that closes a modification annotation (e.g., ']' or ')').</param>
/// <param name="NTermSeparator">Optional separator after N-terminal modification (e.g., "-"). 
/// Null if format doesn't use N-terminal separators.</param>
/// <param name="CTermSeparator">Optional separator before C-terminal modification (e.g., "-").
/// Null if format doesn't use C-terminal separators.</param>
public readonly record struct SequenceFormatSchema(
    string FormatName,
    char ModOpenBracket,
    char ModCloseBracket,
    string? NTermSeparator = null,
    string? CTermSeparator = null)
{
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
