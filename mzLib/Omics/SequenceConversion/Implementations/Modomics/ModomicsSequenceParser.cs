using System.Collections.Immutable;
using System.Text;
using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Parses MODOMICS one-letter RNA modification sequences.
/// </summary>
public sealed class ModomicsSequenceParser : ISequenceParser
{
    public static ModomicsSequenceParser Instance { get; } = new();

    private static readonly HashSet<char> CanonicalResidues = ['A', 'C', 'G', 'U', 'Y'];

    private ModomicsSequenceParser()
    {
    }

    public string FormatName => ModomicsSequenceFormatSchema.Instance.FormatName;

    public SequenceFormatSchema Schema => ModomicsSequenceFormatSchema.Instance;

    public bool CanParse(string input)
    {
        if (string.IsNullOrWhiteSpace(input))
        {
            return false;
        }

        var hasModificationCode = false;
        foreach (var character in input)
        {
            if (char.IsWhiteSpace(character) || CanonicalResidues.Contains(character) || character == 'P')
            {
                continue;
            }

            if (!Mods.ModomicsLoadReport.ModificationsByAbbreviation.ContainsKey(character.ToString()))
            {
                return false;
            }

            hasModificationCode = true;
        }

        return hasModificationCode;
    }

    public CanonicalSequence? Parse(
        string input,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
    {
        warnings ??= new ConversionWarnings();
        if (string.IsNullOrWhiteSpace(input))
        {
            return SequenceConversionHelpers.HandleParserError(
                warnings, mode, ConversionFailureReason.InvalidSequence,
                "Input sequence is null or empty.");
        }

        var baseSequence = new StringBuilder();
        var modifications = new List<CanonicalModification>();

        for (var inputIndex = 0; inputIndex < input.Length; inputIndex++)
        {
            var code = input[inputIndex];
            if (char.IsWhiteSpace(code))
            {
                continue;
            }

            if (CanonicalResidues.Contains(code))
            {
                baseSequence.Append(code);
                continue;
            }

            if (code == 'P')
            {
                baseSequence.Append('Y');
                continue;
            }

            var followingCharacter = inputIndex + 1 < input.Length ? input[inputIndex + 1] : (char?)null;
            var positionType = ModificationPositionType.Residue;
            var residueIndex = baseSequence.Length;

            // A residue code replaces the residue it represents, so the next input
            // character is normally the following residue, not the code's target.
            // Only use the following character to disambiguate multi-target codes.
            if (!ModomicsModificationLookup.Instance.TryResolveCode(code, null, out var modification))
            {
                ModomicsModificationLookup.Instance.TryResolveCode(code, followingCharacter, out modification);
            }

            if (modification is not null)
            {
                var target = modification!.Target?.Motif?.FirstOrDefault() ?? '\0';
                if (IsFivePrimeModification(modification))
                {
                    if (baseSequence.Length != 0)
                    {
                        return HandleCodeError(warnings, mode, code,
                            "A 5' terminal MODOMICS code must occur at the beginning of the sequence.");
                    }

                    positionType = ModificationPositionType.NTerminus;
                    residueIndex = -1;
                }
                else
                {
                    baseSequence.Append(target);
                }

                modifications.Add(new CanonicalModification(
                    positionType,
                    positionType == ModificationPositionType.Residue ? residueIndex : null,
                    target,
                    code.ToString(),
                    modification.MonoisotopicMass,
                    modification.ChemicalFormula,
                    MzLibId: modification.IdWithMotif,
                    MzLibModification: modification));
                continue;
            }

            return HandleCodeError(warnings, mode, code, "The MODOMICS code could not be resolved.");
        }

        if (baseSequence.Length == 0)
        {
            return SequenceConversionHelpers.HandleParserError(
                warnings, mode, ConversionFailureReason.InvalidSequence,
                "No valid sequence characters found.");
        }

        return new CanonicalSequence(
            baseSequence.ToString(),
            modifications.ToImmutableArray(),
            FormatName);
    }

    private static bool IsFivePrimeModification(Modification modification) =>
        modification.LocationRestriction?.Contains("5'-terminal", StringComparison.OrdinalIgnoreCase) == true;

    private static CanonicalSequence? HandleCodeError(
        ConversionWarnings warnings,
        SequenceConversionHandlingMode mode,
        char code,
        string detail)
    {
        warnings.AddIncompatibleItem(code.ToString());
        return SequenceConversionHelpers.HandleParserError(
            warnings, mode, ConversionFailureReason.IncompatibleModifications,
            $"Unable to parse MODOMICS code '{code}': {detail}");
    }
}
