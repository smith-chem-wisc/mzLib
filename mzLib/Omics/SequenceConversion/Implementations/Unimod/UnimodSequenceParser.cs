using System.Globalization;
using System.Text.RegularExpressions;

namespace Omics.SequenceConversion;

public class UnimodSequenceParser : SequenceParserBase
{
    private static readonly Regex LabeledUnimodTokenPattern = new(@"^UNIMOD:\s*\d+$", RegexOptions.Compiled | RegexOptions.IgnoreCase);
    private static readonly Regex UnimodTokenPattern = new(@"^(?:UNIMOD:\s*)?(\d+)$", RegexOptions.Compiled | RegexOptions.IgnoreCase);

    /// <summary>
    /// Singleton instance for convenience.
    /// </summary>
    public static UnimodSequenceParser Instance { get; } = new();

    /// <inheritdoc />
    public override string FormatName => UnimodSequenceFormatSchema.Instance.FormatName;

    /// <inheritdoc />
    public override SequenceFormatSchema Schema => UnimodSequenceFormatSchema.Instance;

    /// <inheritdoc />
    public override bool CanParse(string input)
    {
        if (string.IsNullOrWhiteSpace(input))
            return false;

        if (!input.Contains('[') || !input.Contains(']'))
            return false;

        if (!AreBracketsBalanced(input))
            return false;

        int bracketStart = input.IndexOf('[');
        while (bracketStart >= 0)
        {
            int bracketEnd = FindClosingBracket(input, bracketStart);
            if (bracketEnd == -1)
                return false;

            string modToken = input.Substring(bracketStart + 1, bracketEnd - bracketStart - 1).Trim();
            if (LabeledUnimodTokenPattern.IsMatch(modToken))
                return true;

            bracketStart = input.IndexOf('[', bracketEnd + 1);
        }

        return false;
    }

    /// <summary>
    /// Parses a modification string extracted from brackets.
    /// Expects UNIMOD IDs such as "UNIMOD:35" or "35".
    /// </summary>
    protected override CanonicalModification? ParseModificationString(
        string modString,
        ModificationPositionType positionType,
        int? residueIndex,
        char? targetResidue,
        ConversionWarnings warnings,
        SequenceConversionHandlingMode mode)
    {
        var match = UnimodTokenPattern.Match(modString.Trim());
        if (!match.Success)
        {
            return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                $"Invalid UNIMOD token format: '{modString}'") == null
                ? null
                : default;
        }

        if (!int.TryParse(match.Groups[1].Value, NumberStyles.Integer, CultureInfo.InvariantCulture, out int unimodId))
        {
            return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                $"Failed to parse UNIMOD ID from token: '{modString}'") == null
                ? null
                : default;
        }

        return new CanonicalModification(
            positionType,
            residueIndex,
            targetResidue,
            modString,
            UnimodId: unimodId);
    }
}
