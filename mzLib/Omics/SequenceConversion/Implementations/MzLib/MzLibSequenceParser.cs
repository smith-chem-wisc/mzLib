using System.Text.RegularExpressions;
using Omics.Modifications;


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

        var modification = FindCatalogModification(modString, positionType, residueIndex);

        return new CanonicalModification(
            positionType,
            residueIndex,
            extractedResidue ?? targetResidue,
            modString,
            UnimodId: modification is not null && CanonicalModification.GetUnimodId(modification) is int unimodId
                      && Mods.MatchesUnimodRecordMass(modification, unimodId) ? unimodId : null,
            MzLibId: mzLibId,
            MzLibModification: modification);
    }

    // Entries can share a name: prefer the one of the written type, then the one whose location restriction fits the position.
    private static Modification? FindCatalogModification(string name, ModificationPositionType positionType, int? residueIndex)
    {
        var separator = name.IndexOf(':');
        if (separator <= 0 || separator == name.Length - 1)
            return null;

        var type = name[..separator].Trim();
        var id = name[(separator + 1)..].Trim();

        var candidates = new List<Modification>();
        if (Mods.AllKnownProteinModsDictionary.TryGetValue(id, out var proteinEntry))
            candidates.Add(proteinEntry);
        if (Mods.AllKnownRnaModsDictionary.TryGetValue(id, out var rnaEntry))
            candidates.Add(rnaEntry);
        foreach (var entry in CatalogById.Value[id])
        {
            if (!candidates.Any(c => ReferenceEquals(c, entry)))
                candidates.Add(entry);
        }

        return candidates
            .OrderBy(m => m.ModificationType == type ? 0 : 1)
            .ThenBy(m => PositionRank(m, positionType, residueIndex))
            .FirstOrDefault();
    }

    private static int PositionRank(Modification modification, ModificationPositionType positionType, int? residueIndex)
    {
        var nTerminal = ModificationLocalization.IsNTerminal(modification);
        var cTerminal = ModificationLocalization.IsCTerminal(modification);
        return positionType switch
        {
            ModificationPositionType.NTerminus => nTerminal ? 0 : 1,
            ModificationPositionType.CTerminus => cTerminal ? 0 : 1,
            _ when !nTerminal && !cTerminal => 0,
            _ => (residueIndex == 0 ? nTerminal : cTerminal) ? 1 : 2
        };
    }

    private static readonly Lazy<ILookup<string, Modification>> CatalogById = new(() =>
        Mods.AllProteinModsList.Concat(Mods.AllRnaModsList)
            .Where(m => !string.IsNullOrEmpty(m.IdWithMotif))
            .ToLookup(m => m.IdWithMotif));
}
