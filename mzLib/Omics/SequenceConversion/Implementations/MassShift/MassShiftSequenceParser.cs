using System.Text.RegularExpressions;


namespace Omics.SequenceConversion;

/// <summary>
/// Parses mass shift-format sequence strings into <see cref="CanonicalSequence"/>.
/// 
/// Supports:
/// - Unmodified sequences: "PEPTIDE"
/// - Residue modifications by mass: "PEP[+15.995]TIDE"
/// - Terminal modifications by mass: "[+42.011]PEPTIDE" and "PEPTIDE-[+0.984]"
/// - Negative mass shifts: "PEP[-18.011]TIDE"
/// - Multiple modifications: "PEP[+15.995]TID[+79.966]E"
/// </summary>
public class MassShiftSequenceParser : SequenceParserBase
{
    // Regex to detect mass shift annotations like [+15.995] or [-17.027]
    private static readonly Regex MassShiftPattern = new(@"^[+-]?\d+(\.\d+)?$", RegexOptions.Compiled);

    /// <summary>
    /// Singleton instance for convenience.
    /// </summary>
    public static MassShiftSequenceParser Instance { get; } = new();

    /// <inheritdoc />
    public override string FormatName => MassShiftSequenceFormatSchema.Instance.FormatName;

    /// <inheritdoc />
    public override SequenceFormatSchema Schema => MassShiftSequenceFormatSchema.Instance;

    /// <inheritdoc />
    public override bool CanParse(string input)
    {
        if (string.IsNullOrWhiteSpace(input))
            return false;

        // Mass shift format uses square brackets with numeric content
        if (!input.Contains('[') || !input.Contains(']'))
            return false; // No brackets, not this format

        // Check that brackets are balanced
        if (!AreBracketsBalanced(input))
            return false;

        // Check if any bracket content looks like a mass shift
        int bracketStart = input.IndexOf('[');
        while (bracketStart >= 0)
        {
            int bracketEnd = FindClosingBracket(input, bracketStart);
            if (bracketEnd == -1)
                return false;

            string content = input.Substring(bracketStart + 1, bracketEnd - bracketStart - 1);
            if (MassShiftPattern.IsMatch(content))
                return true; // At least one mass shift found

            bracketStart = input.IndexOf('[', bracketEnd + 1);
        }

        return false; // No mass shifts found
    }

    /// <summary>
    /// Parses a modification string extracted from brackets.
    /// Expects numeric mass shift format (e.g., "+15.995", "-18.011").
    /// </summary>
    protected override CanonicalModification? ParseModificationString(
        string modString,
        ModificationPositionType positionType,
        int? residueIndex,
        char? targetResidue,
        ConversionWarnings warnings,
        SequenceConversionHandlingMode mode)
    {
        // Validate mass shift format
        if (!MassShiftPattern.IsMatch(modString))
        {
            return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                $"Invalid mass shift format: '{modString}'") == null 
                ? null 
                : default;
        }

        // Parse the mass value
        if (!double.TryParse(modString, out double mass))
        {
            return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                $"Failed to parse mass value: '{modString}'") == null 
                ? null 
                : default;
        }

        return new CanonicalModification(
            positionType,
            residueIndex,
            targetResidue,
            modString,
            mass);
    }
}
