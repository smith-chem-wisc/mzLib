using System;
using System.Collections.Generic;
using System.Collections.Concurrent;
using System.Linq;
using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Resolves modifications to their UniProt representations by searching the UniProt PTM database.
/// </summary>
public class UniProtModificationLookup : ModificationLookupBase
{
    public static UniProtModificationLookup Instance { get; } = new();

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

    public UniProtModificationLookup(IEnumerable<Modification>? candidateSet = null, double massTolerance = 0.001)
        : base(candidateSet ?? Mods.UniprotModifications, massTolerance)
    {
    }

    public override string Name => "UniProt";

    public override CanonicalModification? TryResolve(CanonicalModification mod)
    {
        var resolved = base.TryResolve(mod);
        if (resolved == null || string.IsNullOrWhiteSpace(mod.OriginalRepresentation))
        {
            return resolved;
        }

        var normalized = NormalizeRepresentation(mod.OriginalRepresentation);
        var trimmed = TrimUniProtSuffixes(normalized);
        var variants = BuildNameVariants(trimmed, mod.TargetResidue);
        if (variants.Count == 0)
        {
            return resolved;
        }

        var mzLibMod = resolved.Value.MzLibModification;
        if (mzLibMod == null)
        {
            return resolved;
        }

        if (MatchesAnyVariant(mzLibMod, variants))
        {
            return resolved;
        }

        return null;
    }

    public override CanonicalModification? TryResolve(string originalRepresentation, char? targetResidue = null, Chemistry.ChemicalFormula? chemicalFormula = null)
    {
        var resolved = base.TryResolve(originalRepresentation, targetResidue, chemicalFormula);
        if (resolved == null || string.IsNullOrWhiteSpace(originalRepresentation))
        {
            return resolved;
        }

        var normalized = NormalizeRepresentation(originalRepresentation);
        var trimmed = TrimUniProtSuffixes(normalized);
        var variants = BuildNameVariants(trimmed, targetResidue);
        if (variants.Count == 0)
        {
            return resolved;
        }

        var mzLibMod = resolved.Value.MzLibModification;
        if (mzLibMod == null)
        {
            return resolved;
        }

        if (MatchesAnyVariant(mzLibMod, variants))
        {
            return resolved;
        }

        return null;
    }

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

        // If no name matches found, return empty to let other filters (formula, mass) be tried
        // Don't fall back to all candidates - that defeats the purpose of name filtering
        var candidates = nameMatches.Count > 0
            ? nameMatches
            : formulaMatches.Count > 0
                ? formulaMatches
                : Enumerable.Empty<Modification>().ToList();
        
        if (candidates.Count == 0)
        {
            return Enumerable.Empty<Modification>();
        }

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
        
        // Phosphorylation -> Phospho[residue] mapping
        if (trimmedName.IndexOf("phosphory", StringComparison.OrdinalIgnoreCase) >= 0 &&
            ResidueNameMap.TryGetValue(residue, out var residueName))
        {
            var lowerResidueName = residueName.ToLowerInvariant();
            var phosphoName = $"Phospho{lowerResidueName}";
            variants.Add(phosphoName);
            variants.Add($"{phosphoName} on {residue}");
        }
        
        // Oxidation -> specific UniProt names mapping
        if (trimmedName.IndexOf("oxid", StringComparison.OrdinalIgnoreCase) >= 0 &&
            ResidueNameMap.TryGetValue(residue, out var oxidResidueName))
        {
            // Map "Oxidation on X" to UniProt-specific names
            switch (residue)
            {
                case 'M':
                    variants.Add("Methionine sulfoxide");
                    variants.Add("Methionine sulfoxide on M");
                    break;
                case 'W':
                    variants.Add("Oxindolylalanine");
                    variants.Add("Kynurenine");
                    break;
                case 'C':
                    variants.Add("Cysteine sulfenic acid");
                    variants.Add("S-cysteinyl cysteine");
                    break;
                case 'H':
                    variants.Add("2-oxohistidine");
                    break;
                case 'P':
                    variants.Add("Hydroxyproline");
                    break;
                case 'K':
                    variants.Add("Hydroxylysine");
                    break;
                case 'Y':
                    variants.Add("Dihydroxyphenylalanine");
                    break;
            }
        }
        
        // Carbamidomethyl -> S-carbamoylmethylcysteine mapping
        if (trimmedName.IndexOf("carbamidomethyl", StringComparison.OrdinalIgnoreCase) >= 0 && residue == 'C')
        {
            variants.Add("S-carbamoylmethylcysteine");
            variants.Add("S-carbamoylmethylcysteine on C");
        }
        
        // Acetylation -> specific UniProt names mapping
        if (trimmedName.IndexOf("acetyl", StringComparison.OrdinalIgnoreCase) >= 0)
        {
            if (ResidueNameMap.TryGetValue(residue, out var acetylResidueName))
            {
                variants.Add($"N-acetyl{acetylResidueName.ToLowerInvariant()}");
                variants.Add($"N-acetyl{acetylResidueName.ToLowerInvariant()} on {residue}");
            }
            // Also add N6-acetyllysine for internal K acetylation
            if (residue == 'K')
            {
                variants.Add("N6-acetyllysine");
                variants.Add("N6-acetyllysine on K");
            }
            // Also try generic N-terminal acetyl
            variants.Add("N-acetylated residue");
        }
        
        // Malonylation on K -> N6-malonyllysine
        if (trimmedName.IndexOf("malonyl", StringComparison.OrdinalIgnoreCase) >= 0 && residue == 'K')
        {
            variants.Add("N6-malonyllysine");
            variants.Add("N6-malonyllysine on K");
        }
        
        // Hydroxylation -> residue-specific UniProt names
        if (trimmedName.IndexOf("hydroxyl", StringComparison.OrdinalIgnoreCase) >= 0)
        {
            switch (residue)
            {
                case 'K':
                    variants.Add("5-hydroxylysine");
                    variants.Add("5-hydroxylysine on K");
                    break;
                case 'P':
                    variants.Add("Hydroxyproline");
                    variants.Add("Hydroxyproline on P");
                    break;
                case 'N':
                    variants.Add("3-hydroxyasparagine");
                    variants.Add("3-hydroxyasparagine on N");
                    break;
            }
        }
        
        // HexNAc -> O-linked/N-linked (HexNAc) residue names
        if (trimmedName.IndexOf("hexnac", StringComparison.OrdinalIgnoreCase) >= 0)
        {
            if (ResidueNameMap.TryGetValue(residue, out var hexNAcResidueName))
            {
                var lowerName = hexNAcResidueName.ToLowerInvariant();
                switch (residue)
                {
                    case 'T':
                        variants.Add("O-linked (HexNAc) threonine");
                        variants.Add("O-linked (HexNAc) threonine on T");
                        break;
                    case 'S':
                        variants.Add("O-linked (HexNAc) serine");
                        variants.Add("O-linked (HexNAc) serine on S");
                        break;
                    case 'N':
                        variants.Add("N-linked (HexNAc) asparagine");
                        variants.Add("N-linked (HexNAc) asparagine on N");
                        break;
                }
            }
        }
        
        // Methylation -> residue-specific UniProt names
        if (trimmedName.IndexOf("methyl", StringComparison.OrdinalIgnoreCase) >= 0 &&
            trimmedName.IndexOf("dimethyl", StringComparison.OrdinalIgnoreCase) < 0 &&
            trimmedName.IndexOf("trimethyl", StringComparison.OrdinalIgnoreCase) < 0)
        {
            switch (residue)
            {
                case 'K':
                    variants.Add("N6-methyllysine");
                    variants.Add("N6-methyllysine on K");
                    break;
                case 'H':
                    variants.Add("Methylhistidine");
                    variants.Add("Methylhistidine on H");
                    break;
                case 'R':
                    variants.Add("5-methylarginine");
                    variants.Add("5-methylarginine on R");
                    variants.Add("Omega-N-methylarginine");
                    break;
            }
        }
        
        // Dimethylation -> residue-specific UniProt names
        if (trimmedName.IndexOf("dimethyl", StringComparison.OrdinalIgnoreCase) >= 0)
        {
            switch (residue)
            {
                case 'R':
                    variants.Add("Dimethylated arginine");
                    variants.Add("Dimethylated arginine on R");
                    variants.Add("Omega-N,N-dimethylarginine");
                    variants.Add("Asymmetric dimethylarginine");
                    break;
                case 'K':
                    variants.Add("N6,N6-dimethyllysine");
                    variants.Add("N6,N6-dimethyllysine on K");
                    break;
            }
        }
        
        // Trimethylation on K -> N6,N6,N6-trimethyllysine
        if (trimmedName.IndexOf("trimethyl", StringComparison.OrdinalIgnoreCase) >= 0 && residue == 'K')
        {
            variants.Add("N6,N6,N6-trimethyllysine");
            variants.Add("N6,N6,N6-trimethyllysine on K");
        }
        
        // Hydroxybutyrylation on K -> N6-(2-hydroxyisobutyryl)lysine
        if (trimmedName.IndexOf("hydroxybutyryl", StringComparison.OrdinalIgnoreCase) >= 0 && residue == 'K')
        {
            variants.Add("N6-(2-hydroxyisobutyryl)lysine");
            variants.Add("N6-(2-hydroxyisobutyryl)lysine on K");
        }
        
        // Nitrosylation -> residue-specific UniProt names
        if (trimmedName.IndexOf("nitrosyl", StringComparison.OrdinalIgnoreCase) >= 0)
        {
            switch (residue)
            {
                case 'Y':
                    variants.Add("3'-nitrotyrosine");
                    variants.Add("3'-nitrotyrosine on Y");
                    variants.Add("3-nitrotyrosine");
                    break;
                case 'C':
                    variants.Add("S-nitrosocysteine");
                    variants.Add("S-nitrosocysteine on C");
                    break;
            }
        }
        
        // Formylation on K -> N6-formyllysine
        if (trimmedName.IndexOf("formyl", StringComparison.OrdinalIgnoreCase) >= 0 && residue == 'K')
        {
            variants.Add("N6-formyllysine");
            variants.Add("N6-formyllysine on K");
        }
        
        // Crotonylation on K -> N6-crotonyllysine
        if (trimmedName.IndexOf("crotonyl", StringComparison.OrdinalIgnoreCase) >= 0 && residue == 'K')
        {
            variants.Add("N6-crotonyllysine");
            variants.Add("N6-crotonyllysine on K");
        }
        
        // Dehydrobutyrine on T -> 2,3-didehydrobutyrine
        if (trimmedName.IndexOf("dehydrobutyrine", StringComparison.OrdinalIgnoreCase) >= 0 && residue == 'T')
        {
            variants.Add("2,3-didehydrobutyrine");
            variants.Add("2,3-didehydrobutyrine on T");
        }
        
        // Dehydroalanine -> residue-specific UniProt names
        if (trimmedName.IndexOf("dehydroalanine", StringComparison.OrdinalIgnoreCase) >= 0)
        {
            switch (residue)
            {
                case 'S':
                    variants.Add("2,3-didehydroalanine (Ser)");
                    variants.Add("2,3-didehydroalanine (Ser) on S");
                    variants.Add("2, 3-didehydroalanine (Ser)");
                    break;
                case 'C':
                    variants.Add("2,3-didehydroalanine (Cys)");
                    variants.Add("2,3-didehydroalanine (Cys) on C");
                    break;
            }
        }
        
        // Citrullination on R -> Citrulline
        if (trimmedName.IndexOf("citrullin", StringComparison.OrdinalIgnoreCase) >= 0 && residue == 'R')
        {
            variants.Add("Citrulline");
            variants.Add("Citrulline on R");
        }
        
        // Carboxylation -> residue-specific UniProt names
        if (trimmedName.IndexOf("carboxyl", StringComparison.OrdinalIgnoreCase) >= 0)
        {
            switch (residue)
            {
                case 'E':
                    variants.Add("4-carboxyglutamate");
                    variants.Add("4-carboxyglutamate on E");
                    break;
                case 'K':
                    variants.Add("N6-carboxylysine");
                    variants.Add("N6-carboxylysine on K");
                    break;
            }
        }
        
        // Succinylation on K -> N6-succinyllysine
        if (trimmedName.IndexOf("succinyl", StringComparison.OrdinalIgnoreCase) >= 0 && residue == 'K')
        {
            variants.Add("N6-succinyllysine");
            variants.Add("N6-succinyllysine on K");
        }
        
        // Pyridoxal phosphate on K -> N6-(pyridoxal phosphate)lysine
        if (trimmedName.IndexOf("pyridoxal", StringComparison.OrdinalIgnoreCase) >= 0 && residue == 'K')
        {
            variants.Add("N6-(pyridoxal phosphate)lysine");
            variants.Add("N6-(pyridoxal phosphate)lysine on K");
        }
        
        // Sulfonation on Y -> Sulfotyrosine
        if (trimmedName.IndexOf("sulfon", StringComparison.OrdinalIgnoreCase) >= 0 && residue == 'Y')
        {
            variants.Add("Sulfotyrosine");
            variants.Add("Sulfotyrosine on Y");
        }
        
        // Butyrylation on K -> N6-butyryllysine
        if (trimmedName.IndexOf("butyryl", StringComparison.OrdinalIgnoreCase) >= 0 &&
            trimmedName.IndexOf("hydroxybutyryl", StringComparison.OrdinalIgnoreCase) < 0 && residue == 'K')
        {
            variants.Add("N6-butyryllysine");
            variants.Add("N6-butyryllysine on K");
        }
        
        // Glutarylation on K -> N6-glutaryllysine
        if (trimmedName.IndexOf("glutaryl", StringComparison.OrdinalIgnoreCase) >= 0 && residue == 'K')
        {
            variants.Add("N6-glutaryllysine");
            variants.Add("N6-glutaryllysine on K");
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

    private static bool MatchesAnyVariant(Modification mod, IReadOnlyList<string> variants)
    {
        if (mod == null || variants == null || variants.Count == 0)
        {
            return false;
        }

        return MatchesUniProtName(mod.IdWithMotif, variants) || MatchesUniProtName(mod.OriginalId, variants);
    }

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
            .OrderByDescending(mod => HasExactVariantMatch(mod, variants))
            .ThenBy(mod => GetFirstVariantIndex(mod, variants))
            .ThenByDescending(mod => GetBestOverlapScore(mod.IdWithMotif, variants))
            .ThenByDescending(mod => targetResidue != null && mod.Target != null && mod.Target.ToString().Contains(targetResidue.Value))
            .ThenByDescending(mod => mod.ModificationType?.Length ?? 0)
            .ThenBy(mod => mod.IdWithMotif.Length)
            .FirstOrDefault();
    }

    private static bool HasExactVariantMatch(Modification mod, IReadOnlyList<string> variants)
    {
        if (mod == null || variants == null || variants.Count == 0)
        {
            return false;
        }

        foreach (var variant in variants)
        {
            if (string.IsNullOrWhiteSpace(variant))
            {
                continue;
            }

            if (MatchesExactVariant(mod, variant))
            {
                return true;
            }
        }

        return false;
    }

    private static int GetFirstVariantIndex(Modification mod, IReadOnlyList<string> variants)
    {
        if (mod == null || variants == null || variants.Count == 0)
        {
            return int.MaxValue;
        }

        for (int i = 0; i < variants.Count; i++)
        {
            var variant = variants[i];
            if (string.IsNullOrWhiteSpace(variant))
            {
                continue;
            }

            if (MatchesExactVariant(mod, variant) || ContainsVariant(mod.IdWithMotif, new[] { variant }))
            {
                return i;
            }
        }

        return int.MaxValue;
    }

    private static bool MatchesExactVariant(Modification mod, string variant)
    {
        if (mod == null || string.IsNullOrWhiteSpace(variant))
        {
            return false;
        }

        if (!string.IsNullOrEmpty(mod.OriginalId) &&
            string.Equals(mod.OriginalId, variant, StringComparison.OrdinalIgnoreCase))
        {
            return true;
        }

        if (!string.IsNullOrEmpty(mod.IdWithMotif) &&
            string.Equals(mod.IdWithMotif, variant, StringComparison.OrdinalIgnoreCase))
        {
            return true;
        }

        if (!string.IsNullOrEmpty(mod.IdWithMotif) &&
            mod.IdWithMotif.EndsWith(" on X", StringComparison.OrdinalIgnoreCase))
        {
            var baseId = mod.IdWithMotif[..^5].TrimEnd();
            return string.Equals(baseId, variant, StringComparison.OrdinalIgnoreCase);
        }

        return false;
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
