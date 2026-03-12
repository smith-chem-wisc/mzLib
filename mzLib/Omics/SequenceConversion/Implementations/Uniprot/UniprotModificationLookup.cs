using System;
using System.Collections.Generic;
using System.Collections.Concurrent;
using System.Linq;
using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Resolves modifications to their UniProt representations by searching the UniProt PTM database.
/// </summary>
public class UniprotModificationLookup : ModificationLookupBase
{
    public static UniprotModificationLookup Instance { get; } = new();

    private readonly ConcurrentDictionary<(string Name, char? Residue), Modification?> _nameResidueCache = new();

    private static readonly Dictionary<char, string> ResidueNameMap = new()
    {
        ['A'] = "Alanine",
        ['R'] = "Arginine",
        ['N'] = "Asparagine",
        ['D'] = "Aspartic acid",
        ['C'] = "Cysteine",
        ['E'] = "Glutamic acid",
        ['Q'] = "Glutamine",
        ['G'] = "Glycine",
        ['H'] = "Histidine",
        ['I'] = "Isoleucine",
        ['L'] = "Leucine",
        ['K'] = "Lysine",
        ['M'] = "Methionine",
        ['F'] = "Phenylalanine",
        ['P'] = "Proline",
        ['S'] = "Serine",
        ['T'] = "Threonine",
        ['W'] = "Tryptophan",
        ['Y'] = "Tyrosine",
        ['V'] = "Valine",
        ['U'] = "Selenocysteine",
        ['O'] = "Pyrrolysine"
    };

    public UniprotModificationLookup(IEnumerable<Modification>? candidateSet = null, double massTolerance = 0.001)
        : base(candidateSet ?? Mods.UniprotModifications, massTolerance)
    {
    }

    public override string Name => "UniProt";

    protected override IEnumerable<Modification> FilterByName(IEnumerable<Modification> source, string name, char? targetResidue)
    {
        var normalized = NormalizeRepresentation(name);
        if (string.IsNullOrEmpty(normalized))
        {
            return Enumerable.Empty<Modification>();
        }

        var cacheKey = (normalized, targetResidue);
        if (_nameResidueCache.TryGetValue(cacheKey, out var cached) && cached != null)
        {
            return new[] { cached };
        }

        var formulaMatches = new List<Modification>();
        var nameMatches = new List<Modification>();

        source ??= CandidateSet;
        var trimmedName = TrimUniProtSuffixes(normalized);
        var nameVariants = BuildNameVariants(trimmedName, targetResidue);

        foreach (var mod in source)
        {
            if (mod.ChemicalFormula != null && ContainsVariant(mod.IdWithMotif, nameVariants))
            {
                formulaMatches.Add(mod);
            }

            if (MatchesUniProtName(mod.IdWithMotif, nameVariants))
            {
                nameMatches.Add(mod);
            }
        }

        var motifMatches = source.Where(m =>
            m.Target != null &&
            (targetResidue == null ||
             m.Target.ToString().Contains(targetResidue.Value) ||
             m.Target.ToString().Equals("X", StringComparison.OrdinalIgnoreCase)))
            .ToList();

        var candidates = nameMatches.Count > 0
            ? nameMatches
            : formulaMatches.Count > 0
                ? formulaMatches
                : source.ToList();

        var intersected = candidates.Intersect(motifMatches).ToList();
        if (intersected.Count == 0)
        {
            intersected = candidates;
        }

        var resolved = ScoreAndSelect(intersected, nameVariants, targetResidue);
        if (resolved != null)
        {
            _nameResidueCache[cacheKey] = resolved;
            return new[] { resolved };
        }

        return intersected;
    }

    protected override string NormalizeRepresentation(string representation)
    {
        var normalized = base.NormalizeRepresentation(representation);
        if (string.IsNullOrEmpty(normalized))
        {
            return normalized;
        }

        var colonIndex = normalized.IndexOf(':');
        if (colonIndex >= 0)
        {
            var prefix = normalized[..colonIndex].Trim();
            var remainder = normalized[(colonIndex + 1)..].Trim();

            if (prefix.Equals("UniProt", StringComparison.OrdinalIgnoreCase) || prefix.Contains(' '))
            {
                return string.IsNullOrEmpty(remainder) ? normalized : remainder;
            }
        }

        return normalized;
    }

    private static string TrimUniProtSuffixes(string name)
    {
        if (string.IsNullOrEmpty(name))
        {
            return name;
        }

        var trimmed = name.Split(new[] { "-L-" }, StringSplitOptions.None)[0];
        return trimmed;
    }

    private static IReadOnlyList<string> BuildNameVariants(string trimmedName, char? targetResidue)
    {
        var variants = new List<string>();
        if (!string.IsNullOrWhiteSpace(trimmedName))
        {
            variants.Add(trimmedName);
        }

        AddGeneralAtionVariants(trimmedName, variants);
        AddResidueSpecificVariants(trimmedName, targetResidue, variants);

        if (variants.Count == 0)
        {
            return Array.Empty<string>();
        }

        var deduped = new List<string>(variants.Count);
        var seen = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        foreach (var variant in variants)
        {
            if (string.IsNullOrWhiteSpace(variant))
            {
                continue;
            }

            if (seen.Add(variant))
            {
                deduped.Add(variant);
            }
        }

        return deduped;
    }

    private static void AddGeneralAtionVariants(string trimmedName, List<string> variants)
    {
        if (string.IsNullOrEmpty(trimmedName) ||
            trimmedName.IndexOf("ation", StringComparison.OrdinalIgnoreCase) < 0)
        {
            return;
        }

        variants.Add(ReplaceCaseInsensitive(trimmedName, "ation", string.Empty));
        variants.Add(ReplaceCaseInsensitive(trimmedName, "ation", "yl"));
    }

    private static void AddResidueSpecificVariants(string trimmedName, char? targetResidue, List<string> variants)
    {
        if (string.IsNullOrEmpty(trimmedName) || !targetResidue.HasValue)
        {
            return;
        }

        var residue = char.ToUpperInvariant(targetResidue.Value);
        if (trimmedName.IndexOf("phosphory", StringComparison.OrdinalIgnoreCase) >= 0 &&
            ResidueNameMap.TryGetValue(residue, out var residueName))
        {
            var lowerResidueName = residueName.ToLowerInvariant();
            var phosphoName = $"Phospho{lowerResidueName}";
            variants.Add(phosphoName);
            variants.Add($"{phosphoName} on {residue}");
        }
    }

    private static string ReplaceCaseInsensitive(string input, string search, string replacement)
    {
        if (string.IsNullOrEmpty(input) || string.IsNullOrEmpty(search))
        {
            return input;
        }

        var start = input.IndexOf(search, StringComparison.OrdinalIgnoreCase);
        if (start < 0)
        {
            return input;
        }

        return string.Concat(input.AsSpan(0, start), replacement, input.AsSpan(start + search.Length));
    }

    private static bool MatchesUniProtName(string idWithMotif, IReadOnlyList<string> nameVariants) =>
        ContainsVariant(idWithMotif, nameVariants);

    private static bool ContainsVariant(string text, IReadOnlyList<string> variants)
    {
        if (string.IsNullOrEmpty(text) || variants == null || variants.Count == 0)
        {
            return false;
        }

        foreach (var variant in variants)
        {
            if (!string.IsNullOrWhiteSpace(variant) && text.Contains(variant, StringComparison.OrdinalIgnoreCase))
            {
                return true;
            }
        }

        return false;
    }

    private Modification? ScoreAndSelect(IList<Modification> candidates, IReadOnlyList<string> nameVariants, char? targetResidue)
    {
        if (candidates.Count == 0)
        {
            return null;
        }

        var nonLabel = candidates.Where(c => !c.IdWithMotif.StartsWith("Label", StringComparison.OrdinalIgnoreCase)).ToList();
        if (nonLabel.Count > 0)
        {
            candidates = nonLabel;
        }

        var variants = nameVariants?.Count > 0 ? nameVariants : Array.Empty<string>();

        return candidates
            .OrderByDescending(mod => GetBestOverlapScore(mod.IdWithMotif, variants))
            .ThenByDescending(mod => targetResidue != null && mod.Target != null && mod.Target.ToString().Contains(targetResidue.Value))
            .ThenByDescending(mod => mod.ModificationType?.Length ?? 0)
            .ThenBy(mod => mod.IdWithMotif.Length)
            .FirstOrDefault();
    }

    private static int GetBestOverlapScore(string idWithMotif, IReadOnlyList<string> variants)
    {
        if (string.IsNullOrEmpty(idWithMotif) || variants == null || variants.Count == 0)
        {
            return 0;
        }

        var best = 0;
        foreach (var variant in variants)
        {
            if (string.IsNullOrWhiteSpace(variant))
            {
                continue;
            }

            best = Math.Max(best, GetOverlapScore(idWithMotif, variant));
        }

        return best;
    }
}
