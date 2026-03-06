using System.Collections.Immutable;
using System.Text;
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
public class MassShiftSequenceParser : ISequenceParser
{
    // Regex to detect mass shift annotations like [+15.995] or [-17.027]
    private static readonly Regex MassShiftPattern = new(@"^[+-]?\d+(\.\d+)?$", RegexOptions.Compiled);

    /// <summary>
    /// Singleton instance for convenience.
    /// </summary>
    public static MassShiftSequenceParser Instance { get; } = new();

    /// <inheritdoc />
    public string FormatName => MassShiftSequenceFormatSchema.Instance.FormatName;

    /// <inheritdoc />
    public SequenceFormatSchema Schema => MassShiftSequenceFormatSchema.Instance;

    /// <inheritdoc />
    public bool CanParse(string input)
    {
        if (string.IsNullOrWhiteSpace(input))
            return false;

        // Mass shift format uses square brackets with numeric content
        if (!input.Contains('[') || !input.Contains(']'))
            return false; // No brackets, not this format

        // Check that brackets are balanced
        int depth = 0;
        bool hasSquareBrackets = false;
        foreach (char c in input)
        {
            if (c == '[')
            {
                depth++;
                hasSquareBrackets = true;
            }
            else if (c == ']')
            {
                depth--;
            }

            if (depth < 0)
                return false; // More closing than opening
        }

        if (!hasSquareBrackets || depth != 0)
            return false; // Unbalanced brackets

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

        // Check for N-terminal modification: "[+mass]sequence"
        if (input.StartsWith("["))
        {
            int closeBracket = FindClosingBracket(input, 0);
            if (closeBracket == -1)
            {
                return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                    "Unclosed bracket at start of sequence.");
            }

            string modString = input.Substring(1, closeBracket - 1);
            if (!MassShiftPattern.IsMatch(modString))
            {
                return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                    $"Invalid mass shift format at N-terminus: '{modString}'");
            }

            double mass = double.Parse(modString);
            var nTermMod = new CanonicalModification(
                ModificationPositionType.NTerminus,
                null,
                null,
                modString,
                mass);
            modifications.Add(nTermMod);

            // Skip past the closing bracket (no separator for N-term)
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
                if (!MassShiftPattern.IsMatch(modString))
                {
                    return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                        $"Invalid mass shift format at position {i}: '{modString}'");
                }

                double mass = double.Parse(modString);
                char? targetResidue = residueIndex > 0 ? baseSequence[residueIndex - 1] : null;

                var mod = new CanonicalModification(
                    ModificationPositionType.Residue,
                    residueIndex - 1,
                    targetResidue,
                    modString,
                    mass);
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
                    if (!MassShiftPattern.IsMatch(modString))
                    {
                        return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                            $"Invalid mass shift format at C-terminus: '{modString}'");
                    }

                    double mass = double.Parse(modString);
                    char? targetResidue = baseSequence.Length > 0 ? baseSequence[baseSequence.Length - 1] : null;

                    var cTermMod = new CanonicalModification(
                        ModificationPositionType.CTerminus,
                        null,
                        targetResidue,
                        modString,
                        mass);
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
            MassShiftSequenceFormatSchema.Instance.FormatName);
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
