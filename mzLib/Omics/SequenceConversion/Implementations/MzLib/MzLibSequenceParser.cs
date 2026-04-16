using System.Text.RegularExpressions;


namespace Omics.SequenceConversion;

/// <summary>
/// Parses mzLib-format sequence strings into <see cref="CanonicalSequence"/>.
/// 
/// Supports:
/// - Unmodified sequences: "PEPTIDE"
/// - Residue modifications: "PEP[Oxidation on M]TIDE"
/// - Terminal modifications: "[Acetyl]PEPTIDE" and "PEPTIDE-[Amidated]"
/// - Multiple modifications: "PEP[Oxidation on M]TID[Phospho on S]E"
/// 
/// Note: For mass shift notation (e.g., "[+15.995]"), use MassShiftSequenceParser instead.
/// </summary>
public class MzLibSequenceParser : SequenceParserBase
{
    /// <summary>
    /// Singleton instance for convenience.
    /// </summary>
    public static MzLibSequenceParser Instance { get; } = new();

    /// <inheritdoc />
    public override string FormatName => MzLibSequenceFormatSchema.Instance.FormatName;

    /// <inheritdoc />
    public override SequenceFormatSchema Schema => MzLibSequenceFormatSchema.Instance;

    /// <inheritdoc />
    public override bool CanParse(string input)
    {
        if (string.IsNullOrWhiteSpace(input))
            return false;

        if (input.Contains("unimod:", StringComparison.InvariantCultureIgnoreCase))
            return false; // unimod: identifiers are not mzLib format

        // mzLib format uses square brackets
        // Check for balanced brackets and valid structure
        bool hasSquareBrackets = input.Contains('[') && input.Contains(']');
        bool hasParentheses = input.Contains('(') && input.Contains(')');

        // If it has parentheses but no square brackets, it's probably not mzLib format
        if (hasParentheses && !hasSquareBrackets)
            return false;

        // Check that brackets are balanced
        return AreBracketsBalanced(input);
    }

    /// <summary>
    /// Parses a modification string extracted from brackets.
    /// Expects mzLib-style modification identifiers (e.g., "Oxidation on M", "Common Fixed:Carbamidomethyl on C").
    /// </summary>
    protected override CanonicalModification? ParseModificationString(
        string modString,
        ModificationPositionType positionType,
        int? residueIndex,
        char? targetResidue,
        ConversionWarnings warnings,
        SequenceConversionHandlingMode mode)
    {
        // mzLib format uses modification identifiers, not mass shifts
        // Mass shifts should be parsed by MassShiftSequenceParser instead
        
        // Check if it has a modification type prefix (e.g., "Common Fixed:Carbamidomethyl on C")
        string mzLibId = modString;
        if (modString.Contains(':'))
        {
            // Keep the full string as-is - the lookup will handle it
            mzLibId = modString;
        }

        // Try to extract target residue from "on X" suffix if not already known
        char? extractedResidue = targetResidue;
        if (!extractedResidue.HasValue)
        {
            var onMatch = Regex.Match(modString, @"\s+on\s+(\w)$", RegexOptions.IgnoreCase);
            if (onMatch.Success)
            {
                extractedResidue = onMatch.Groups[1].Value[0];
            }
        }

        return new CanonicalModification(
            positionType,
            residueIndex,
            extractedResidue ?? targetResidue,
            modString,
            MzLibId: mzLibId);
    }
}
