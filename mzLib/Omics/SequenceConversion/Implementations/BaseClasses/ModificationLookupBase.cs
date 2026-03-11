using Chemistry;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace Omics.SequenceConversion;

/// <summary>
/// Base class for modification lookups providing shared resolution logic.
/// Subclasses implement database-specific lookup strategies.
/// </summary>
public abstract class ModificationLookupBase : IModificationLookup
{
    protected IReadOnlyCollection<Modification> CandidateSet { get; }
    protected double? MassTolerance { get; }

    protected ModificationLookupBase(IEnumerable<Modification>? candidateSet, double? massTolerance)
    {
        IReadOnlyCollection<Modification>? resolvedCandidates = candidateSet as IReadOnlyCollection<Modification>;
        if (resolvedCandidates == null && candidateSet != null)
        {
            resolvedCandidates = candidateSet.ToList();
        }

        CandidateSet = resolvedCandidates ?? Array.Empty<Modification>();
        MassTolerance = massTolerance;
    }

    /// <inheritdoc />
    public abstract string Name { get; }

    /// <inheritdoc />
    public CanonicalModification? TryResolve(CanonicalModification mod)
    {
        // If already resolved, return as-is
        if (mod.IsResolved)
            return mod;

        // Try primary resolution strategy (database-specific identifiers)
        var primary = ResolveFromCandidates(GetPrimaryCandidates(mod), mod);
        if (primary != null)
        {
            return mod.WithResolvedModification(primary, mod.ResidueIndex, mod.PositionType);
        }

        // Try to resolve by original representation
        var byName = ResolveFromCandidates(FilterByName(CandidateSet, mod.OriginalRepresentation, mod.TargetResidue), mod);
        if (byName != null)
        {
            return mod.WithResolvedModification(byName, mod.ResidueIndex, mod.PositionType);
        }

        // Try chemical formula fallback
        if (mod.ChemicalFormula != null)
        {
            var formulaMatch = ResolveFromCandidates(FilterByFormula(CandidateSet, mod.ChemicalFormula), mod);
            if (formulaMatch != null)
            {
                return mod.WithResolvedModification(formulaMatch, mod.ResidueIndex, mod.PositionType);
            }
        }

        // Try mass-based fallback (optional, requires tolerance)
        if (mod.MonoisotopicMass.HasValue)
        {
            var massMatch = ResolveFromCandidates(FilterByMass(CandidateSet, mod.MonoisotopicMass.Value), mod);
            if (massMatch != null)
            {
                return mod.WithResolvedModification(massMatch, mod.ResidueIndex, mod.PositionType);
            }
        }

        return null;
    }

    /// <inheritdoc />
    public CanonicalModification? TryResolve(string originalRepresentation, char? targetResidue = null, ChemicalFormula? chemicalFormula = null)
    {
        var candidates = FilterByName(CandidateSet, originalRepresentation, targetResidue);
        var resolved = SelectUniqueUsingEvidence(candidates, targetResidue, chemicalFormula, null);
        if (resolved != null)
        {
            return CreateCanonicalModification(originalRepresentation, resolved, targetResidue);
        }

        if (chemicalFormula != null)
        {
            var formulaCandidates = FilterByFormula(CandidateSet, chemicalFormula);
            resolved = SelectUniqueUsingEvidence(formulaCandidates, targetResidue, chemicalFormula, null);
            if (resolved != null)
            {
                return CreateCanonicalModification(originalRepresentation, resolved, targetResidue);
            }
        }

        return null;
    }

    /// <summary>
    /// Primary resolution strategy using database-specific identifiers.
    /// Implementations should return the set of potential matches; the base class
    /// will apply additional filters to disambiguate.
    /// </summary>
    protected virtual IEnumerable<Modification> GetPrimaryCandidates(CanonicalModification mod) => Enumerable.Empty<Modification>();

    /// <summary>
    /// Applies normalization/expansion to a name and filters the provided candidates.
    /// </summary>
    protected virtual IEnumerable<Modification> FilterByName(IEnumerable<Modification> source, string name, char? targetResidue)
    {
        if (string.IsNullOrWhiteSpace(name))
            return Enumerable.Empty<Modification>();

        var normalized = NormalizeRepresentation(name);
        if (string.IsNullOrEmpty(normalized))
            return Enumerable.Empty<Modification>();

        source ??= CandidateSet;

        HashSet<string> potentialStrings = new(StringComparer.OrdinalIgnoreCase) { normalized };
        foreach (var candidate in ExpandNameCandidates(normalized, targetResidue))
        {
            potentialStrings.Add(candidate);
        }

        return source.Where(modification => potentialStrings.Any(candidate => MatchesIdentifier(modification, candidate)));
    }

    /// <summary>
    /// Filters by chemical formula.
    /// </summary>
    protected virtual IEnumerable<Modification> FilterByFormula(IEnumerable<Modification> source, ChemicalFormula formula) =>
        (source ?? CandidateSet).Where(m => m.ChemicalFormula != null && m.ChemicalFormula.Equals(formula));

    /// <summary>
    /// Filters by mass within the configured tolerance.
    /// </summary>
    protected virtual IEnumerable<Modification> FilterByMass(IEnumerable<Modification> source, double mass)
    {
        if (!MassTolerance.HasValue)
            return Enumerable.Empty<Modification>();

        var tolerance = MassTolerance.Value;
        return (source ?? CandidateSet).Where(m => m.MonoisotopicMass.HasValue &&
                                                   Math.Abs(m.MonoisotopicMass.Value - mass) <= tolerance);
    }

    /// <summary>
    /// Filters by identifiers that match the standard IdWithMotif/OriginalId/Type:Id patterns.
    /// </summary>
    protected virtual IEnumerable<Modification> FilterByIdentifier(IEnumerable<Modification> source, string identifier)
    {
        if (string.IsNullOrWhiteSpace(identifier))
            return Enumerable.Empty<Modification>();

        var normalized = NormalizeRepresentation(identifier);
        if (string.IsNullOrEmpty(normalized))
            return Enumerable.Empty<Modification>();

        source ??= CandidateSet;
        return source.Where(m => MatchesIdentifier(m, normalized));
    }

    /// <summary>
    /// Filters candidates by UNIMOD database references.
    /// </summary>
    protected IEnumerable<Modification> FilterByUnimodId(IEnumerable<Modification> source, int unimodId)
    {
        source ??= CandidateSet;
        var idString = unimodId.ToString(CultureInfo.InvariantCulture);

        return source.Where(m =>
            (m.DatabaseReference != null &&
             m.DatabaseReference.Any(kvp => kvp.Key.Equals("UNIMOD", StringComparison.OrdinalIgnoreCase) &&
                                             kvp.Value.Any(value => value.Contains(idString, StringComparison.OrdinalIgnoreCase))))
            || (!string.IsNullOrEmpty(m.Accession) && m.Accession.Contains(idString, StringComparison.OrdinalIgnoreCase)));
    }

    /// <summary>
    /// Creates a CanonicalModification from a resolved Modification.
    /// Uses WithResolvedModification to extract UNIMOD ID and other metadata.
    /// </summary>
    protected CanonicalModification CreateCanonicalModification(
        string originalRepresentation,
        Modification resolvedMod,
        char? targetResidue)
    {
        // Create a minimal unresolved CanonicalModification, then resolve it
        var unresolved = new CanonicalModification(
            PositionType: ModificationPositionType.Residue, // Default, caller should adjust if terminal
            ResidueIndex: null,
            TargetResidue: targetResidue,
            OriginalRepresentation: originalRepresentation);

        // WithResolvedModification extracts UNIMOD ID and other metadata from the resolved mod
        return unresolved.WithResolvedModification(resolvedMod);
    }

    #region Candidate resolution helpers

    private Modification? ResolveFromCandidates(IEnumerable<Modification> candidates, CanonicalModification evidence)
    {
        return SelectUniqueUsingEvidence(candidates, evidence.TargetResidue, evidence.ChemicalFormula, evidence.MonoisotopicMass);
    }

    protected Modification? SelectUniqueUsingEvidence(
        IEnumerable<Modification> candidates,
        char? targetResidue,
        ChemicalFormula? formula,
        double? mass)
    {
        var current = (candidates as IList<Modification>) ?? candidates?.ToList() ?? new List<Modification>();

        if (current.Count == 0)
            return null;

        var unique = SelectUnique(current, targetResidue);
        if (unique != null)
            return unique;

        if (formula != null)
        {
            current = current.Where(m => m.ChemicalFormula != null && m.ChemicalFormula.Equals(formula)).ToList();
            unique = SelectUnique(current, targetResidue);
            if (unique != null)
                return unique;
        }

        if (mass.HasValue)
        {
            var massFiltered = FilterByMass(current, mass.Value).ToList();
            unique = SelectUnique(massFiltered, targetResidue);
            if (unique != null)
                return unique;
        }

        return null;
    }

    private static Modification? SelectUnique(IEnumerable<Modification> candidates, char? targetResidue)
    {
        var list = (candidates as IList<Modification>) ?? candidates.ToList();
        if (list.Count == 1)
        {
            return list[0];
        }

        if (list.Count > 1 && targetResidue.HasValue)
        {
            var residueMatches = list.Where(m =>
                m.Target != null &&
                m.Target.ToString().Contains(targetResidue.Value)).ToList();

            if (residueMatches.Count == 1)
            {
                return residueMatches[0];
            }
        }

        return null;
    }

    /// <summary>
    /// Calculates the overlap score between the modification ID with motif and the trimmed name.
    /// The score represents the length of the longest common substring between the two strings.
    /// </summary>
    /// <param name="idWithMotif">The modification ID with motif.</param>
    /// <param name="trimmedName">The trimmed name of the modification.</param>
    /// <returns>The overlap score, which is the length of the longest common substring.</returns>
    protected static int GetOverlapScore(string idWithMotif, string trimmedName)
    {
        int overlapScore = 0;
        for (int i = 0; i < idWithMotif.Length; i++)
        {
            for (int j = 0; j < trimmedName.Length; j++)
            {
                int k = 0;
                while (i + k < idWithMotif.Length && j + k < trimmedName.Length && idWithMotif[i + k] == trimmedName[j + k])
                {
                    k++;
                }
                overlapScore = Math.Max(overlapScore, k);
            }
        }
        return overlapScore;
    }

    #endregion

    #region Name Normalization and Expansion 
    protected virtual string NormalizeRepresentation(string representation) => representation?.Trim() ?? string.Empty;

    protected virtual IEnumerable<string> ExpandNameCandidates(string normalizedRepresentation, char? targetResidue)
    {
        if (string.IsNullOrEmpty(normalizedRepresentation))
            yield break;

        yield return normalizedRepresentation;
        if (targetResidue.HasValue)
            yield return $"{normalizedRepresentation} on {targetResidue.Value}";

        bool containsNTerminalRepresentation = normalizedRepresentation.Contains("on N-terminus", StringComparison.OrdinalIgnoreCase);
        if (containsNTerminalRepresentation)
            foreach (var additionalName in ExpandNTerminus(normalizedRepresentation, targetResidue))
                yield return additionalName;

        bool isIsobaric = IsIsobaric(normalizedRepresentation);
        if (isIsobaric)
            foreach (var additionalName in ExpandIsobaricTags(normalizedRepresentation, targetResidue))
                yield return additionalName;


        if (!normalizedRepresentation.Contains(":")) 
            yield break;

        var second = normalizedRepresentation.Split(':')[1];
        yield return second;
        if (targetResidue.HasValue)
            yield return $"{second} on {targetResidue.Value}";

        if (containsNTerminalRepresentation)
            foreach (var additionalName in ExpandNTerminus(second, targetResidue))
                yield return additionalName;
        if (isIsobaric)
            foreach (var additionalName in ExpandIsobaricTags(second, targetResidue))
                yield return additionalName;
    }

    private IEnumerable<string> ExpandNTerminus(string input, char? targetResidue)
    {
        yield return input.Replace("on N-terminus", "on X", StringComparison.OrdinalIgnoreCase);
        if (targetResidue.HasValue)
            yield return input.Replace("on N-terminus", $"on {targetResidue.Value}", StringComparison.OrdinalIgnoreCase);
    }

    private bool IsIsobaric(string input) => input.Contains("TMT", StringComparison.OrdinalIgnoreCase)
        || input.Contains("iTRAQ", StringComparison.OrdinalIgnoreCase)
        || input.Contains("plex", StringComparison.OrdinalIgnoreCase)
        || input.Contains("DiLeu", StringComparison.OrdinalIgnoreCase);

    private IEnumerable<string> ExpandIsobaricTags(string input, char? targetResidue)
    {
        // Handle Isobaric Labeling variations (e.g., TMT, iTRAQ) that may be represented with or without "-plex".
        if (input.Contains("-plex", StringComparison.OrdinalIgnoreCase))
        {
            yield return input.Replace("-plex", "plex", StringComparison.OrdinalIgnoreCase);
        }
        else if (input.Contains("plex", StringComparison.OrdinalIgnoreCase))
        {
            yield return input.Replace("plex", "-plex", StringComparison.OrdinalIgnoreCase);
        }

        if (input.Contains("TMT18", StringComparison.OrdinalIgnoreCase))
        {
            yield return input.Replace("TMT18", "TMTpro", StringComparison.OrdinalIgnoreCase);
        }
        else if (input.Contains("TMTpro", StringComparison.OrdinalIgnoreCase))
        {
            yield return input.Replace("TMTpro", "TMT18", StringComparison.OrdinalIgnoreCase);
        }

        if (input.Contains("TMT10pro", StringComparison.OrdinalIgnoreCase))
        {
            yield return input.Replace("TMT10pro", "TMT10", StringComparison.OrdinalIgnoreCase);
        }
        else if (input.Contains("TMT10", StringComparison.OrdinalIgnoreCase))
        {
            yield return input.Replace("TMT10", "TMT10pro", StringComparison.OrdinalIgnoreCase);
        }
    }

    #endregion

    protected virtual bool MatchesIdentifier(Modification modification, string identifier) =>
        modification.IdWithMotif == identifier ||
        modification.OriginalId == identifier ||
        $"{modification.ModificationType}:{modification.IdWithMotif}" == identifier;

}
