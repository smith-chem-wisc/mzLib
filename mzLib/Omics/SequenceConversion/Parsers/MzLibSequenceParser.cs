using System.Collections.Immutable;
using System.Text;
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
public class MzLibSequenceParser : ISequenceParser
{

    /// <summary>
    /// Singleton instance for convenience.
    /// </summary>
    public static MzLibSequenceParser Instance { get; } = new();

    /// <inheritdoc />
    public string FormatName => MzLibSequenceFormatSchema.Instance.FormatName;

    /// <inheritdoc />
    public SequenceFormatSchema Schema => MzLibSequenceFormatSchema.Instance;

    /// <inheritdoc />
    public bool CanParse(string input)
    {
        if (string.IsNullOrWhiteSpace(input))
            return false;

        // mzLib format uses square brackets
        // Check for balanced brackets and valid structure
        bool hasSquareBrackets = input.Contains('[') && input.Contains(']');
        bool hasParentheses = input.Contains('(') && input.Contains(')');

        // If it has parentheses but no square brackets, it's probably not mzLib format
        if (hasParentheses && !hasSquareBrackets)
            return false;

        // Check that brackets are balanced
        int depth = 0;
        foreach (char c in input)
        {
            if (c == '[') depth++;
            else if (c == ']') depth--;

            if (depth < 0)
                return false; // More closing than opening
        }

        return depth == 0;
    }

    /// <inheritdoc />
    public CanonicalSequence? Parse(
        string input,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
    {
        warnings ??= new ConversionWarnings();

        if (string.IsNullOrWhiteSpace(input))
        {
            return HandleError(warnings, mode, ConversionFailureReason.InvalidSequence,
                "Input sequence is null or empty.");
        }

        try
        {
            return ParseInternal(input, warnings, mode);
        }
        catch (SequenceConversionException)
        {
            throw; // Re-throw our own exceptions
        }
        catch (Exception ex)
        {
            return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                $"Unexpected error parsing sequence: {ex.Message}");
        }
    }

    private CanonicalSequence? ParseInternal(string input, ConversionWarnings warnings, SequenceConversionHandlingMode mode)
    {
        var baseSequence = new StringBuilder();
        var modifications = new List<CanonicalModification>();
        int residueIndex = 0;
        int i = 0;

        // Check for N-terminal modification: "[mod]sequence" (no separator)
        if (input.StartsWith("["))
        {
            int closeBracket = FindClosingBracket(input, 0);
            if (closeBracket == -1)
            {
                return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                    "Unclosed bracket at start of sequence.");
            }

            string modString = input.Substring(1, closeBracket - 1);
            var nTermMod = ParseModificationString(modString, ModificationPositionType.NTerminus, null, null);
            modifications.Add(nTermMod);

            // Skip past the closing bracket (no separator for N-term in mzLib format)
            i = closeBracket + 1;
        }

        // Parse the main sequence
        while (i < input.Length)
        {
            char c = input[i];

            if (c == '[')
            {
                // This is a modification on the previous residue
                int closeBracket = FindClosingBracket(input, i);
                if (closeBracket == -1)
                {
                    return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                        $"Unclosed bracket at position {i}.");
                }

                string modString = input.Substring(i + 1, closeBracket - i - 1);
                char? targetResidue = residueIndex > 0 ? baseSequence[residueIndex - 1] : null;

                var mod = ParseModificationString(modString, ModificationPositionType.Residue, residueIndex - 1, targetResidue);
                modifications.Add(mod);

                i = closeBracket + 1;
            }
            else if (c == '-')
            {
                // Check if this is a C-terminal modification separator
                if (i + 1 < input.Length && input[i + 1] == '[')
                {
                    int closeBracket = FindClosingBracket(input, i + 1);
                    if (closeBracket == -1)
                    {
                        return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                            $"Unclosed bracket for C-terminal modification at position {i + 1}.");
                    }

                    string modString = input.Substring(i + 2, closeBracket - i - 2);
                    char? targetResidue = baseSequence.Length > 0 ? baseSequence[baseSequence.Length - 1] : null;

                    var cTermMod = ParseModificationString(modString, ModificationPositionType.CTerminus, null, targetResidue);
                    modifications.Add(cTermMod);

                    i = closeBracket + 1;
                }
                else
                {
                    // Unexpected hyphen in sequence - treat as error
                    return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                        $"Unexpected hyphen at position {i}.");
                }
            }
            else if (char.IsLetter(c))
            {
                // Regular amino acid/nucleotide
                baseSequence.Append(c);
                residueIndex++;
                i++;
            }
            else if (char.IsWhiteSpace(c))
            {
                // Skip whitespace
                i++;
            }
            else
            {
                // Unknown character
                warnings.AddWarning($"Unexpected character '{c}' at position {i}, skipping.");
                i++;
            }
        }

        if (baseSequence.Length == 0)
        {
            return HandleError(warnings, mode, ConversionFailureReason.InvalidSequence,
                "No valid sequence characters found.");
        }

        return new CanonicalSequence(
            baseSequence.ToString(),
            modifications.ToImmutableArray(),
            MzLibSequenceFormatSchema.Instance.FormatName);
    }

    /// <summary>
    /// Parses a modification string extracted from brackets.
    /// Expects mzLib-style modification identifiers (e.g., "Oxidation on M", "Common Fixed:Carbamidomethyl on C").
    /// </summary>
    private CanonicalModification ParseModificationString(
        string modString,
        ModificationPositionType positionType,
        int? residueIndex,
        char? targetResidue)
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

    /// <summary>
    /// Finds the matching closing bracket for an opening bracket.
    /// </summary>
    private static int FindClosingBracket(string input, int openPos)
    {
        if (input[openPos] != '[')
            return -1;

        int depth = 1;
        for (int i = openPos + 1; i < input.Length; i++)
        {
            if (input[i] == '[') depth++;
            else if (input[i] == ']')
            {
                depth--;
                if (depth == 0)
                    return i;
            }
        }
        return -1;
    }

    /// <summary>
    /// Handles an error based on the handling mode.
    /// </summary>
    private static CanonicalSequence? HandleError(
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
            SequenceConversionHandlingMode.ReturnNull => null,
            _ => null
        };
    }
}
