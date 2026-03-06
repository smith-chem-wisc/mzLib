using System.Text.RegularExpressions;

namespace Omics.SequenceConversion;

/// <summary>
/// Applies regex replacements sequentially to produce the target representation.
/// </summary>
public abstract class RegexBasedConverter : SequenceConverterBase
{
    protected abstract IReadOnlyList<(Regex Pattern, string Replacement)> Patterns { get; }

    protected override string ConvertSequenceCore(string fullSequence, SequenceConversionHandlingMode handlingMode,
        IList<string> warnings)
    {
        string converted = fullSequence;
        foreach (var (pattern, replacement) in Patterns)
        {
            converted = pattern.Replace(converted, replacement);
        }

        return converted;
    }
}
