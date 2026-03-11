using System;
using System.Collections.Generic;
using System.Linq;
using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Resolves modifications to their UniProt representations by searching the UniProt PTM database.
/// </summary>
public class UniprotModificationLookup : ModificationLookupBase
{
    public static UniprotModificationLookup Instance { get; } = new();

    public UniprotModificationLookup(IEnumerable<Modification>? candidateSet = null, double massTolerance = 0.001)
        : base(candidateSet ?? Mods.UniprotModifications, massTolerance)
    {
    }

    public override string Name => "UniProt";

    protected override IEnumerable<Modification> GetPrimaryCandidates(CanonicalModification mod)
    {
        if (mod.MzLibModification == null)
        {
            return Enumerable.Empty<Modification>();
        }

        if (IsUniProt(mod.MzLibModification))
        {
            return new[] { mod.MzLibModification };
        }

        var residue = mod.TargetResidue ?? mod.MzLibModification.Target?.ToString().FirstOrDefault();
        return FindCandidatesFromSource(mod.MzLibModification, residue);
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

    private static bool IsUniProt(Modification modification) =>
        modification.ModificationType != null &&
        modification.ModificationType.Equals("UniProt", StringComparison.OrdinalIgnoreCase);

    private IEnumerable<Modification> FindCandidatesFromSource(Modification source, char? residue)
    {
        var results = new List<Modification>();

        void AddMatches(IEnumerable<Modification> matches)
        {
            foreach (var match in matches)
            {
                if (!results.Contains(match))
                {
                    results.Add(match);
                }
            }
        }

        if (source.ChemicalFormula != null)
        {
            AddMatches(FilterByFormula(CandidateSet, source.ChemicalFormula));
        }

        if (source.MonoisotopicMass.HasValue)
        {
            AddMatches(FilterByMass(CandidateSet, source.MonoisotopicMass.Value));
        }

        foreach (var identifier in EnumerateIdentifiers(source))
        {
            AddMatches(FilterByIdentifier(CandidateSet, identifier));
        }

        return results;
    }

    private static IEnumerable<string> EnumerateIdentifiers(Modification source)
    {
        if (!string.IsNullOrWhiteSpace(source.IdWithMotif))
        {
            yield return source.IdWithMotif;
        }

        if (!string.IsNullOrWhiteSpace(source.OriginalId))
        {
            yield return source.OriginalId;
        }

        if (!string.IsNullOrWhiteSpace(source.ModificationType) && !string.IsNullOrWhiteSpace(source.IdWithMotif))
        {
            yield return $"{source.ModificationType}:{source.IdWithMotif}";
        }
    }
}
