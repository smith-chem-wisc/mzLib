using Chemistry;
using Omics.Modifications;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace Omics.SequenceConversion;

/// <summary>
/// Base class for modification lookups providing shared resolution logic.
/// Subclasses implement database-specific lookup strategies.
/// </summary>
public abstract class ModificationLookupBase : IModificationLookup
{
    #region State and Construction

    private readonly ConcurrentDictionary<string, Modification?> _representationCache = new(StringComparer.OrdinalIgnoreCase);

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

    #endregion

    #region Public API

    /// <inheritdoc />
    public abstract string Name { get; }

    /// <inheritdoc />
    public CanonicalModification? TryResolve(CanonicalModification mod)
    {
        if (mod.IsResolved)
        {
            return mod;
        }

        var normalized = NormalizeRepresentation(mod.OriginalRepresentation);
        var resolved = ResolveWithCache(
            normalized,
            mod.TargetResidue,
            mod.ChemicalFormula,
            mod.MonoisotopicMass,
            context => ResolveInternal(context, () => GetPrimaryCandidates(mod)));

        return resolved != null
            ? mod.WithResolvedModification(resolved, mod.ResidueIndex, mod.PositionType)
            : null;
    }

    /// <inheritdoc />
    public CanonicalModification? TryResolve(string originalRepresentation, char? targetResidue = null, ChemicalFormula? chemicalFormula = null)
    {
        var normalized = NormalizeRepresentation(originalRepresentation);
        var resolved = ResolveWithCache(
            normalized,
            targetResidue,
            chemicalFormula,
            mass: null,
            context => ResolveInternal(context, null));

        return resolved != null
            ? CreateCanonicalModification(originalRepresentation, resolved, targetResidue)
            : null;
    }

    #endregion

    #region Caching

    private Modification? ResolveWithCache(
        string normalizedRepresentation,
        char? residue,
        ChemicalFormula? formula,
        double? mass,
        Func<ResolutionContext, Modification?> resolver)
    {
        var cacheKey = BuildCacheKey(normalizedRepresentation, residue, formula, mass);
        if (TryGetCachedResolution(cacheKey, out var cached))
        {
            return cached;
        }

        var context = new ResolutionContext(normalizedRepresentation, residue, formula, mass);
        var resolved = resolver(context);
        StoreCachedResolution(cacheKey, resolved);
        return resolved;
    }

    private protected bool TryGetCachedResolution(string cacheKey, out Modification? resolved) =>
        _representationCache.TryGetValue(cacheKey, out resolved);

    private protected void StoreCachedResolution(string cacheKey, Modification? resolved) =>
        _representationCache[cacheKey] = resolved;

    private protected string BuildCacheKey(string representation, char? residue, ChemicalFormula? formula, double? mass)
    {
        var builder = new StringBuilder(representation ?? string.Empty);

        if (residue.HasValue)
        {
            builder.Append("|res:");
            builder.Append(char.ToUpperInvariant(residue.Value));
        }

        if (formula != null)
        {
            builder.Append("|cf:");
            builder.Append(formula.Formula);
        }

        if (mass.HasValue)
        {
            builder.Append("|mass:");
            builder.Append(mass.Value.ToString("G6", CultureInfo.InvariantCulture));

            if (MassTolerance.HasValue)
            {
                builder.Append("@tol:");
                builder.Append(MassTolerance.Value.ToString("G6", CultureInfo.InvariantCulture));
            }
        }

        return builder.ToString();
    }

    private readonly struct ResolutionContext
    {
        public ResolutionContext(string representation, char? residue, ChemicalFormula? formula, double? mass)
        {
            Representation = representation ?? string.Empty;
            TargetResidue = residue;
            ChemicalFormula = formula;
            Mass = mass;
        }

        public string Representation { get; }
        public char? TargetResidue { get; }
        public ChemicalFormula? ChemicalFormula { get; }
        public double? Mass { get; }
    }

    #endregion

    #region Resolution Pipeline

    /// <summary>
    /// Primary resolution strategy using database-specific identifiers.
    /// Implementations should return the set of potential matches; the base class
    /// will apply additional filters to disambiguate.
    /// </summary>
    protected virtual IEnumerable<Modification> GetPrimaryCandidates(CanonicalModification mod) => Enumerable.Empty<Modification>();

    private Modification? ResolveInternal(ResolutionContext context, Func<IEnumerable<Modification>>? primaryCandidatesFactory)
    {
        foreach (var stage in EnumerateResolutionStages(context, primaryCandidatesFactory))
        {
            var resolved = SelectUniqueUsingContext(stage, context);
            if (resolved != null)
            {
                return resolved;
            }
        }

        return null;
    }

    private IEnumerable<IEnumerable<Modification>> EnumerateResolutionStages(
        ResolutionContext context,
        Func<IEnumerable<Modification>>? primaryCandidatesFactory)
    {
        if (primaryCandidatesFactory != null)
        {
            var primary = primaryCandidatesFactory();
            if (primary != null)
            {
                yield return primary;
            }
        }

        yield return FilterByName(CandidateSet, context.Representation, context.TargetResidue);

        if (context.ChemicalFormula != null)
        {
            yield return FilterByFormula(CandidateSet, context.ChemicalFormula);
        }

        if (context.Mass.HasValue)
        {
            yield return FilterByMass(CandidateSet, context.Mass.Value);
        }
    }

    private Modification? SelectUniqueUsingContext(IEnumerable<Modification> candidates, ResolutionContext context)
    {
        if (candidates == null)
        {
            return null;
        }

        var working = (candidates as IList<Modification>) ?? candidates.ToList();
        if (working.Count == 0)
        {
            return null;
        }

        return SelectBestCandidate(working, context.TargetResidue);
    }

    private static Modification? SelectBestCandidate(IList<Modification> candidates, char? targetResidue)
    {
        if (candidates.Count == 0)
        {
            return null;
        }

        if (targetResidue.HasValue)
        {
            var residueMatches = candidates
                .Where(m => m.Target != null && m.Target.ToString().IndexOf(targetResidue.Value) >= 0)
                .ToList();

            if (residueMatches.Count > 0)
            {
                return residueMatches[0];
            }
        }

        return candidates.Count == 1 ? candidates[0] : null;
    }

    #endregion

    #region Candidate Filters

    /// <summary>
    /// Applies normalization/expansion to a name and filters the provided candidates.
    /// </summary>
    protected virtual IEnumerable<Modification> FilterByName(IEnumerable<Modification> source, string name, char? targetResidue)
    {
        if (string.IsNullOrWhiteSpace(name))
        {
            return Enumerable.Empty<Modification>();
        }

        var normalized = NormalizeRepresentation(name);
        if (string.IsNullOrEmpty(normalized))
        {
            return Enumerable.Empty<Modification>();
        }

        source ??= CandidateSet;

        var potentialStrings = new HashSet<string>(StringComparer.OrdinalIgnoreCase) { normalized };
        foreach (var candidate in ExpandNameCandidates(normalized, targetResidue))
        {
            potentialStrings.Add(candidate);
        }

        return FilterByIdentifierSet(source, potentialStrings);
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
        {
            return Enumerable.Empty<Modification>();
        }

        var tolerance = MassTolerance.Value;
        return (source ?? CandidateSet).Where(m =>
            m.MonoisotopicMass.HasValue &&
            Math.Abs(m.MonoisotopicMass.Value - mass) <= tolerance);
    }

    /// <summary>
    /// Filters by identifiers that match the standard IdWithMotif/OriginalId/Type:Id patterns.
    /// </summary>
    protected virtual IEnumerable<Modification> FilterByIdentifier(IEnumerable<Modification> source, string identifier)
    {
        if (string.IsNullOrWhiteSpace(identifier))
        {
            return Enumerable.Empty<Modification>();
        }

        var normalized = NormalizeRepresentation(identifier);
        if (string.IsNullOrEmpty(normalized))
        {
            return Enumerable.Empty<Modification>();
        }

        var identifiers = new HashSet<string>(StringComparer.OrdinalIgnoreCase) { normalized };
        return FilterByIdentifierSet(source, identifiers);
    }

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

    private IEnumerable<Modification> FilterByIdentifierSet(IEnumerable<Modification> source, HashSet<string> identifiers)
    {
        if (identifiers == null || identifiers.Count == 0)
        {
            return Enumerable.Empty<Modification>();
        }

        source ??= CandidateSet;
        return source.Where(modification => MatchesIdentifierSet(modification, identifiers));
    }

    private bool MatchesIdentifierSet(Modification modification, HashSet<string> identifiers)
    {
        if (modification == null)
        {
            return false;
        }

        foreach (var identifier in identifiers)
        {
            if (MatchesIdentifier(modification, identifier))
            {
                return true;
            }
        }

        return false;
    }

    #endregion

    #region Canonical Helpers

    /// <summary>
    /// Creates a CanonicalModification from a resolved Modification.
    /// Uses WithResolvedModification to extract UNIMOD ID and other metadata.
    /// </summary>
    protected CanonicalModification CreateCanonicalModification(
        string originalRepresentation,
        Modification resolvedMod,
        char? targetResidue)
    {
        var unresolved = new CanonicalModification(
            PositionType: ModificationPositionType.Residue,
            ResidueIndex: null,
            TargetResidue: targetResidue,
            OriginalRepresentation: originalRepresentation);

        return unresolved.WithResolvedModification(resolvedMod);
    }

    #endregion

    #region Name Normalization and Expansion

    protected virtual string NormalizeRepresentation(string representation) => representation?.Trim() ?? string.Empty;

    protected virtual IEnumerable<string> ExpandNameCandidates(string normalizedRepresentation, char? targetResidue)
    {
        if (string.IsNullOrEmpty(normalizedRepresentation))
        {
            yield break;
        }

        yield return normalizedRepresentation;
        if (targetResidue.HasValue)
        {
            yield return $"{normalizedRepresentation} on {targetResidue.Value}";
        }

        bool containsNTerminalRepresentation = normalizedRepresentation.Contains("on N-terminus", StringComparison.OrdinalIgnoreCase);
        if (containsNTerminalRepresentation)
        {
            foreach (var additionalName in ExpandNTerminus(normalizedRepresentation, targetResidue))
            {
                yield return additionalName;
            }
        }

        bool isIsobaric = IsIsobaric(normalizedRepresentation);
        if (isIsobaric)
        {
            foreach (var additionalName in ExpandIsobaricTags(normalizedRepresentation, targetResidue))
            {
                yield return additionalName;
            }
        }

        var delimiterIndex = normalizedRepresentation.IndexOf(':');
        if (delimiterIndex < 0)
        {
            yield break;
        }

        var suffix = normalizedRepresentation[(delimiterIndex + 1)..].Trim();
        if (suffix.Length == 0)
        {
            yield break;
        }

        yield return suffix;
        if (targetResidue.HasValue)
        {
            yield return $"{suffix} on {targetResidue.Value}";
        }

        if (containsNTerminalRepresentation)
        {
            foreach (var additionalName in ExpandNTerminus(suffix, targetResidue))
            {
                yield return additionalName;
            }
        }

        if (isIsobaric)
        {
            foreach (var additionalName in ExpandIsobaricTags(suffix, targetResidue))
            {
                yield return additionalName;
            }
        }
    }

    private IEnumerable<string> ExpandNTerminus(string input, char? targetResidue)
    {
        yield return input.Replace("on N-terminus", "on X", StringComparison.OrdinalIgnoreCase);
        if (targetResidue.HasValue)
        {
            yield return input.Replace("on N-terminus", $"on {targetResidue.Value}", StringComparison.OrdinalIgnoreCase);
        }
    }

    private bool IsIsobaric(string input) => input.Contains("TMT", StringComparison.OrdinalIgnoreCase)
        || input.Contains("iTRAQ", StringComparison.OrdinalIgnoreCase)
        || input.Contains("plex", StringComparison.OrdinalIgnoreCase)
        || input.Contains("DiLeu", StringComparison.OrdinalIgnoreCase);

    private IEnumerable<string> ExpandIsobaricTags(string input, char? targetResidue)
    {
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

    #region Utilities

    protected virtual bool MatchesIdentifier(Modification modification, string identifier) =>
        modification.IdWithMotif == identifier ||
        modification.OriginalId == identifier ||
        (modification.ModificationType != null && modification.IdWithMotif != null &&
         $"{modification.ModificationType}:{modification.IdWithMotif}" == identifier);

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
}
