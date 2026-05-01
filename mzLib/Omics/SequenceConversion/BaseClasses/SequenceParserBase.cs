using System.Collections.Immutable;
using System.Text;

namespace Omics.SequenceConversion;

/// <summary>
/// Base class for sequence parsers that use bracket notation for modifications.
/// Provides common parsing logic for formats like mzLib and mass-shift.
/// </summary>
public abstract class SequenceParserBase : ISequenceParser
{
    /// <inheritdoc />
    public abstract string FormatName { get; }

    /// <inheritdoc />
    public abstract SequenceFormatSchema Schema { get; }

    /// <summary>
    /// Determines whether a C-terminal separator that is not followed by a modification bracket should be treated as an error.
    /// </summary>
    protected virtual bool ThrowOnDanglingCTermSeparator => true;

    /// <inheritdoc />
    public abstract bool CanParse(string input);

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

    /// <summary>
    /// Internal parsing implementation shared across bracket-based parsers.
    /// </summary>
    protected CanonicalSequence? ParseInternal(string input, ConversionWarnings warnings, SequenceConversionHandlingMode mode)
    {
        var baseSequence = new StringBuilder();
        var modifications = new List<CanonicalModification>();
        int residueIndex = 0;
        int i = 0;

        // Check for N-terminal modification: "[mod]sequence" or "[mod]-sequence"
        if (input.StartsWith(Schema.ModOpenBracket))
        {
            int closeBracket = FindClosingBracket(input, 0);
            if (closeBracket == -1)
            {
                return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                    "Unclosed bracket at start of sequence.");
            }

            string modString = input.Substring(1, closeBracket - 1);
            var nTermMod = ParseModificationString(modString, ModificationPositionType.NTerminus, null, null, warnings, mode);
            if (nTermMod == null)
                return null; // Error already handled

            modifications.Add(nTermMod.Value);

            // Skip past the closing bracket
            i = closeBracket + 1;

            // Check for N-terminal separator (if schema defines one)
            if (Schema.NTermSeparator != null && !string.IsNullOrEmpty(Schema.NTermSeparator))
            {
                if (input.Substring(i).StartsWith(Schema.NTermSeparator))
                {
                    i += Schema.NTermSeparator.Length;
                }
                else
                {
                    return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                        $"Expected N-terminal separator '{Schema.NTermSeparator}' after N-terminal modification at position {i}.");
                }
            }
        }

        // Parse the main sequence
        while (i < input.Length)
        {
            char c = input[i];

            if (c == Schema.ModOpenBracket)
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

                var mod = ParseModificationString(modString, ModificationPositionType.Residue, residueIndex - 1, targetResidue, warnings, mode);
                if (mod == null)
                    return null; // Error already handled

                modifications.Add(mod.Value);

                i = closeBracket + 1;
            }
            else if (Schema.CTermSeparator != null && 
                     !string.IsNullOrEmpty(Schema.CTermSeparator) &&
                     input.Substring(i).StartsWith(Schema.CTermSeparator))
            {
                // Check if this is a C-terminal modification separator
                int separatorEnd = i + Schema.CTermSeparator.Length;
                if (separatorEnd < input.Length && input[separatorEnd] == Schema.ModOpenBracket)
                {
                    int closeBracket = FindClosingBracket(input, separatorEnd);
                    if (closeBracket == -1)
                    {
                        return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                            $"Unclosed bracket for C-terminal modification at position {separatorEnd}.");
                    }

                    string modString = input.Substring(separatorEnd + 1, closeBracket - separatorEnd - 1);
                    char? targetResidue = baseSequence.Length > 0 ? baseSequence[baseSequence.Length - 1] : null;

                    var cTermMod = ParseModificationString(modString, ModificationPositionType.CTerminus, null, targetResidue, warnings, mode);
                    if (cTermMod == null)
                        return null; // Error already handled

                    modifications.Add(cTermMod.Value);

                    i = closeBracket + 1;
                }
                else
                {
                    if (ThrowOnDanglingCTermSeparator)
                    {
                        // Separator found but no bracket follows - treat as error
                        return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                            $"C-terminal separator '{Schema.CTermSeparator}' at position {i} not followed by modification bracket.");
                    }

                    // Legacy-compatible behavior: tolerate separator noise when not part of a terminal modification token.
                    i += Schema.CTermSeparator.Length;
                }
            }
            else if (c == Schema.ModCloseBracket)
            {
                // Unmatched closing bracket
                return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                    $"Unmatched closing bracket at position {i}.");
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
            FormatName);
    }

    /// <summary>
    /// Parses a modification string extracted from brackets.
    /// Must be implemented by derived classes to handle format-specific modification syntax.
    /// </summary>
    /// <param name="modString">The modification string (without brackets).</param>
    /// <param name="positionType">The position type (N-terminus, residue, or C-terminus).</param>
    /// <param name="residueIndex">The zero-based residue index (null for terminal modifications).</param>
    /// <param name="targetResidue">The target residue character (may be null).</param>
    /// <param name="warnings">Warnings collection to add any warnings.</param>
    /// <param name="mode">Handling mode for errors.</param>
    /// <returns>The parsed CanonicalModification, or null if parsing failed.</returns>
    protected abstract CanonicalModification? ParseModificationString(
        string modString,
        ModificationPositionType positionType,
        int? residueIndex,
        char? targetResidue,
        ConversionWarnings warnings,
        SequenceConversionHandlingMode mode);

    /// <summary>
    /// Finds the matching closing bracket for an opening bracket.
    /// </summary>
    protected int FindClosingBracket(string input, int openPos)
    {
        if (input[openPos] != Schema.ModOpenBracket)
            return -1;

        int depth = 1;
        for (int i = openPos + 1; i < input.Length; i++)
        {
            if (input[i] == Schema.ModOpenBracket) depth++;
            else if (input[i] == Schema.ModCloseBracket)
            {
                depth--;
                if (depth == 0)
                    return i;
            }
        }
        return -1;
    }

    /// <summary>
    /// Checks if brackets are balanced in the input string.
    /// </summary>
    protected bool AreBracketsBalanced(string input)
    {
        int depth = 0;
        foreach (char c in input)
        {
            if (c == Schema.ModOpenBracket) depth++;
            else if (c == Schema.ModCloseBracket) depth--;

            if (depth < 0)
                return false; // More closing than opening
        }

        return depth == 0;
    }

    /// <summary>
    /// Handles an error based on the handling mode.
    /// </summary>
    protected static CanonicalSequence? HandleError(
        ConversionWarnings warnings,
        SequenceConversionHandlingMode mode,
        ConversionFailureReason reason,
        string message)
    {
        return SequenceConversionHelpers.HandleParserError(warnings, mode, reason, message);
    }
}
