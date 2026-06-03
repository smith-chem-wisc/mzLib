using Chemistry;
using MzLibUtil;
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

    /// <summary>
    /// Instance-scoped cache to prevent cross-lookup pollution.
    /// Each lookup type maintains its own cache keyed by representation + context.
    /// </summary>
    private readonly ConcurrentDictionary<string, Modification?> _instanceCache = new(StringComparer.OrdinalIgnoreCase);

    protected IReadOnlyCollection<Modification> CandidateSet { get; }
    protected Tolerance MassTolerance { get; }

    protected ModificationLookupBase(IEnumerable<Modification>? candidateSet, double? massTolerance)
    {
        IReadOnlyCollection<Modification>? resolvedCandidates = candidateSet as IReadOnlyCollection<Modification>;
        if (resolvedCandidates == null && candidateSet != null)
        {
            resolvedCandidates = candidateSet.ToList();
        }

        CandidateSet = resolvedCandidates ?? Mods.AllKnownMods;
        MassTolerance = new AbsoluteTolerance(massTolerance ?? 0.01);
    }

    #endregion

    #region Public API

    /// <inheritdoc />
    public abstract string Name { get; }

    /// <inheritdoc />
    public virtual CanonicalModification? TryResolve(CanonicalModification mod)
    {
        // Only skip lookup if the mod is already resolved to THIS lookup's database
        // If it's resolved to a different database (e.g., MetaMorpheus), we still need to find the equivalent
        if (mod.IsResolved && IsResolvedToThisDatabase(mod))
        {
            return mod;
        }

        var normalized = NormalizeRepresentation(mod.OriginalRepresentation);
        var resolved = ResolveWithCache(
            normalized,
            mod.TargetResidue,
            mod.ChemicalFormula,
            mod.MonoisotopicMass ?? mod.ChemicalFormula?.MonoisotopicMass ?? null,
            MassTolerance,
            mod.PositionType,
            context => ResolveInternal(context, () => GetPrimaryCandidates(mod)));

        if (resolved != null)
            return mod.WithResolvedModification(resolved, mod.ResidueIndex, mod.PositionType);

        return null;
    }

    /// <inheritdoc />
    public virtual CanonicalModification? TryResolve(string originalRepresentation, char? targetResidue = null, ChemicalFormula? chemicalFormula = null, ModificationPositionType? positionType = null)
    {
        var normalized = NormalizeRepresentation(originalRepresentation);
        var resolved = ResolveWithCache(
            normalized,
            targetResidue,
            chemicalFormula,
            chemicalFormula?.MonoisotopicMass ?? null,
            MassTolerance,
            positionType, 
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
        Tolerance massTolerance,
        ModificationPositionType? term,
        Func<ResolutionContext, Modification?> resolver)
    {
        var cacheKey = BuildCacheKey(normalizedRepresentation, residue, formula, mass, massTolerance, term);
        if (_instanceCache.TryGetValue(cacheKey, out var cached))
        {
            return cached;
        }

        var context = new ResolutionContext(normalizedRepresentation, residue, formula, mass, term);
        var resolved = resolver(context);
        _instanceCache[cacheKey] = resolved;
        return resolved;
    }

    protected static string BuildCacheKey(string representation, char? residue, ChemicalFormula? formula, double? mass, Tolerance massTolerance, ModificationPositionType? term)
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


            builder.Append("@tol:");
            builder.Append(massTolerance.Value.ToString("G6", CultureInfo.InvariantCulture));
        }

        if (term.HasValue)
        {
            builder.Append("|term:");
            builder.Append(term);
        }

        return builder.ToString();
    }

    private readonly struct ResolutionContext
    {
        public ResolutionContext(string representation, char? residue, ChemicalFormula? formula, double? mass, ModificationPositionType? term)
        {
            Representation = representation ?? string.Empty;
            TargetResidue = residue;
            ChemicalFormula = formula;
            Mass = mass;
            Term = term;
        }

        public string Representation { get; }
        public char? TargetResidue { get; }
        public ChemicalFormula? ChemicalFormula { get; }
        public double? Mass { get; }
        public ModificationPositionType? Term { get; }
    }

    #endregion

    #region Resolution Pipeline

    /// <summary>
    /// Primary resolution strategy using database-specific identifiers.
    /// Implementations should return the set of potential matches; the base class
    /// will apply additional filters to disambiguate.
    /// Returns empty if no specific identifiers (MzLibId/UnimodId) are available.
    /// </summary>
    protected virtual IEnumerable<Modification> GetPrimaryCandidates(CanonicalModification mod)
    {
        if (!string.IsNullOrEmpty(mod.MzLibId))
        {
            return FilterByIdentifier(CandidateSet, mod.MzLibId);
        }

        if (mod.UnimodId.HasValue)
        {
            return FilterByUnimodId(CandidateSet, mod.UnimodId.Value);
        }

        // No specific identifiers - return empty to signal that we should use cumulative filters
        return [];
    }

    /// <summary>
    /// Resolves a modification by applying filters cumulatively.
    /// Strategy:
    /// 1. If primary candidates (from MzLibId/UnimodId) exist, try to select best match
    /// 2. Otherwise, apply cumulative filters (name, formula, mass) to narrow down CandidateSet
    /// 3. If filters don't match anything, return null (don't return arbitrary candidates)
    /// </summary>
    private Modification? ResolveInternal(ResolutionContext context, Func<IEnumerable<Modification>>? primaryCandidatesFactory)
    {
        // Step 1: Check for primary candidates (from specific identifiers like MzLibId/UnimodId)
        if (primaryCandidatesFactory != null)
        {
            var primary = primaryCandidatesFactory();
            var primaryList = primary?.ToList() ?? [];
            if (primaryList.Count > 0)
            {
                // We have primary candidates from specific identifiers - try to select best match
                var uniqueFromPrimary = SelectBestCandidate(primaryList, context.TargetResidue);
                if (uniqueFromPrimary != null)
                {
                    return uniqueFromPrimary;
                }
                // If SelectBestCandidate returned null (ambiguous), continue with cumulative filters
                // using the primary candidates as the starting set
                return ApplyCumulativeFilters(primaryList, context);
            }
        }

        // Step 2: No primary candidates - apply cumulative filters starting from full candidate set
        return ApplyCumulativeFilters(CandidateSet, context);
    }

    /// <summary>
    /// Applies all available filters cumulatively to narrow down candidates.
    /// Returns a unique match if found, or null if no unique match can be determined.
    /// Key behavior: if identifying information is provided (name, formula, mass) but
    /// doesn't match anything, returns null rather than falling through to return
    /// an arbitrary candidate.
    /// </summary>
    private Modification? ApplyCumulativeFilters(IEnumerable<Modification> startingCandidates, ResolutionContext context)
    {
        var candidates = startingCandidates.ToList();
        if (candidates.Count == 0)
        {
            return null;
        }

        // Track whether any identifying filter was applied and matched
        bool anyFilterApplied = false;
        bool anyFilterMatched = false;

        var hasRepresentation = !string.IsNullOrWhiteSpace(context.Representation);
        var isMassRepresentation = hasRepresentation && IsMassRepresentation(context.Representation);

        // Apply name filter if we have a non-mass representation
        if (hasRepresentation && !isMassRepresentation)
        {
            anyFilterApplied = true;
            var byName = FilterByName(candidates, context.Representation, context.TargetResidue).ToList();
            if (byName.Count > 0)
            {
                anyFilterMatched = true;
                candidates = byName;
                var result = SelectBestCandidate(candidates, context.TargetResidue);
                if (result != null)
                {
                    return result;
                }
            }
            else
            {
                if (context.ChemicalFormula == null && !context.Mass.HasValue)
                {
                    return null;
                }
            }
        }

        // Apply formula filter if available
        if (context.ChemicalFormula != null)
        {
            anyFilterApplied = true;
            var byFormula = FilterByFormula(candidates, context.ChemicalFormula).ToList();
            if (byFormula.Count > 0)
            {
                anyFilterMatched = true;
                candidates = byFormula;
                var result = SelectBestCandidate(candidates, context.TargetResidue);
                if (result != null)
                {
                    return result;
                }
            }
            else
            {
                return null;
            }
        }

        // Apply mass filter if available
        if (context.Mass.HasValue)
        {
            anyFilterApplied = true;
            var byMass = FilterByMass(candidates, context.Mass.Value).ToList();
            if (byMass.Count > 0)
            {
                anyFilterMatched = true;
                candidates = byMass;
                var result = SelectBestCandidate(candidates, context.TargetResidue);
                if (result != null)
                {
                    return result;
                }
            }
            else
            {
                return null;
            }
        }

        if (context.TargetResidue.HasValue)
        {
            anyFilterApplied = true;
            var byMotif = FilterByMotif(candidates, context.TargetResidue.Value).ToList();
            if (byMotif.Count > 0)
            {
                anyFilterMatched = true;
                candidates = byMotif;
                var result = SelectBestCandidate(candidates, context.TargetResidue);
                if (result != null)
                {
                    return result;
                }
            }
        }

        // If we applied filters but none matched, return null
        // (don't return an arbitrary candidate when identifying info didn't match)
        if (anyFilterApplied && !anyFilterMatched)
        {
            return null;
        }

        if (context.Term.HasValue && candidates.Count > 0)
        {
            candidates = FilterByTerm(candidates, context.Term.Value).ToList();
        }

        // Final attempt with remaining candidates (only if filters matched or no filters applied)
        return anyFilterMatched ? SelectBestCandidate(candidates, context.TargetResidue) : null;
    }

    /// <summary>
    /// Selects the best candidate from a list based on residue matching and other criteria.
    /// Does NOT mutate the input list.
    /// </summary>
    private static Modification? SelectBestCandidate(IList<Modification> candidates, char? targetResidue)
    {
        if (candidates.Count == 0)
        {
            return null;
        }

        if (candidates.Count == 1)
        {
            return candidates[0];
        }

        // Create a working copy to avoid mutating input
        var working = candidates.ToList();
        var target = targetResidue?.ToString();

        // Filter by residue match if target residue is specified
        if (targetResidue.HasValue && !string.IsNullOrEmpty(target))
        {
            var residueMatches = working.Where(c => 
                c.Target?.Motif == target || 
                c.Target?.Motif == "X" ||  // X is wildcard
                string.IsNullOrEmpty(c.Target?.Motif)).ToList();
            
            if (residueMatches.Count > 0)
            {
                // Prefer exact residue match over wildcard
                var exactMatches = residueMatches.Where(c => c.Target?.Motif == target).ToList();
                if (exactMatches.Count > 0)
                {
                    working = exactMatches;
                }
                else
                {
                    working = residueMatches;
                }
            }
        }

        // Filter out isobaric mods without diagnostic ions
        var nonIsobaricOrWithDiagnostic = working.Where(c =>
            !IsIsobaric(c.OriginalId) || (c.DiagnosticIons != null && c.DiagnosticIons.Count >= 1)).ToList();
        
        if (nonIsobaricOrWithDiagnostic.Count > 0)
        {
            working = nonIsobaricOrWithDiagnostic;
        }

        if (working.Count == 1)
        {
            return working[0];
        }

        // Check if remaining candidates are functionally equivalent (same formula)
        var withFormula = working.Where(m => m.ChemicalFormula != null).ToList();
        var distinctFormulas = withFormula.Select(m => m.ChemicalFormula?.Formula).Distinct().ToList();
        
        if (distinctFormulas.Count == 1 && withFormula.Count > 0)
        {
            // All candidates have the same formula - they're chemically equivalent
            // Prefer the one with the shortest name (generic form over stereospecific)
            // e.g., prefer "Methionine sulfoxide" over "Methionine (R)-sulfoxide"
            return withFormula.OrderBy(m => m.OriginalId?.Length ?? int.MaxValue)
                .ThenBy(m => m.OriginalId) // Stable sort for determinism
                .First();
        }

        // If we still have multiple candidates, prefer those with formulas
        if (withFormula.Count == 1)
        {
            return withFormula[0];
        }

        // Check if all remaining candidates have the same target residue
        // If so, they're functionally equivalent and we can return any
        var distinctByResidue = working.DistinctBy(m => m.Target?.Motif ?? "").ToList();
        if (distinctByResidue.Count == 1)
        {
            return distinctByResidue[0];
        }

        // Multiple distinct candidates with different target residues - truly ambiguous
        // Return null to signal we can't resolve without more context
        return null;
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
            return [];
        }

        var normalized = NormalizeRepresentation(name);
        if (string.IsNullOrEmpty(normalized))
        {
            return [];
        }

        source ??= CandidateSet;

        var potentialStrings = new HashSet<string>(StringComparer.OrdinalIgnoreCase) { normalized };
        foreach (var candidate in ExpandNameCandidates(normalized, targetResidue))
            potentialStrings.Add(candidate);

        int colonIndex = normalized.IndexOf(":", StringComparison.Ordinal);
        if (colonIndex > 0)
        {
            var suffix = normalized[(colonIndex + 1)..];
            foreach (var candidate in ExpandNameCandidates(suffix, targetResidue))
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
        return (source ?? CandidateSet).Where(m => MassTolerance.Within(m.MonoisotopicMass!.Value, mass));
    }

    /// <summary>
    /// Filters by identifiers that match the standard IdWithMotif/OriginalId/Type:Id patterns.
    /// </summary>
    protected virtual IEnumerable<Modification> FilterByIdentifier(IEnumerable<Modification> source, string identifier)
    {
        if (string.IsNullOrWhiteSpace(identifier))
        {
            return [];
        }

        var normalized = NormalizeRepresentation(identifier);
        if (string.IsNullOrEmpty(normalized))
        {
            return [];
        }

        var identifiers = new HashSet<string>(StringComparer.OrdinalIgnoreCase) { normalized };
        return FilterByIdentifierSet(source, identifiers);
    }

    protected IEnumerable<Modification> FilterByMotif(IEnumerable<Modification> source, char motif)
    {
        source ??= CandidateSet;
        return source.Where(m => m.Target != null && 
                                 (m.Target.Motif == motif.ToString() || m.Target.Motif == "X" || string.IsNullOrEmpty(m.Target.Motif)));
    }

    protected virtual IEnumerable<Modification> FilterByTerm(IEnumerable<Modification> source, ModificationPositionType term)
    {
        source ??= CandidateSet;
        return source.Where(m => MatchesTermRestriction(m.LocationRestriction, term));
    }

    protected static bool MatchesTermRestriction(string? locationRestriction, ModificationPositionType term)
    {
        if (string.IsNullOrWhiteSpace(locationRestriction) ||
            locationRestriction.Equals("Anywhere.", StringComparison.OrdinalIgnoreCase) ||
            locationRestriction.Equals("Unassigned.", StringComparison.OrdinalIgnoreCase))
        {
            return true;
        }

        return term switch
        {
            ModificationPositionType.NTerminus => IsNTerminalRestriction(locationRestriction),
            ModificationPositionType.CTerminus => IsCTerminalRestriction(locationRestriction),
            ModificationPositionType.Residue => !IsNTerminalRestriction(locationRestriction) && !IsCTerminalRestriction(locationRestriction),
            _ => true
        };
    }

    private static bool IsNTerminalRestriction(string locationRestriction)
    {
        return locationRestriction.Contains("N-terminal", StringComparison.OrdinalIgnoreCase) ||
               locationRestriction.Contains("5'-terminal", StringComparison.OrdinalIgnoreCase);
    }

    private static bool IsCTerminalRestriction(string locationRestriction)
    {
        return locationRestriction.Contains("C-terminal", StringComparison.OrdinalIgnoreCase) ||
               locationRestriction.Contains("3'-terminal", StringComparison.OrdinalIgnoreCase);
    }

    protected static bool IsMassRepresentation(string representation)
    {
        if (string.IsNullOrWhiteSpace(representation))
        {
            return false;
        }

        var trimmed = representation.Trim().Trim('[', ']');
        trimmed = trimmed.TrimStart('+');
        return double.TryParse(trimmed, NumberStyles.Float, CultureInfo.InvariantCulture, out _);
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
            return [];
        }

        source ??= CandidateSet;
        return source.Where(modification => MatchesIdentifierSet(modification, identifiers))
            .Distinct();
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
                if (containsNTerminalRepresentation)
                    foreach (var isobarExpanded in ExpandNTerminus(additionalName, targetResidue))
                        yield return isobarExpanded;
            }
        }
    }

    private static IEnumerable<string> ExpandNTerminus(string input, char? targetResidue)
    {
        yield return input.Replace("on N-terminus", "on X", StringComparison.OrdinalIgnoreCase);
        if (targetResidue.HasValue)
        {
            yield return input.Replace("on N-terminus", $"on {targetResidue.Value}", StringComparison.OrdinalIgnoreCase);
        }
    }

    private static bool IsIsobaric(string input) => input.Contains("TMT", StringComparison.OrdinalIgnoreCase)
                                                   || input.Contains("iTRAQ", StringComparison.OrdinalIgnoreCase)
                                                   || input.Contains("plex", StringComparison.OrdinalIgnoreCase)
                                                   || input.Contains("DiLeu", StringComparison.OrdinalIgnoreCase);

    private static IEnumerable<string> ExpandIsobaricTags(string input, char? targetResidue)
    {
        if (input.Contains("-plex", StringComparison.OrdinalIgnoreCase))
        {
            yield return input.Replace("-plex", "plex", StringComparison.OrdinalIgnoreCase);
            yield return input.Replace("-plex", "", StringComparison.OrdinalIgnoreCase);
        }
        else if (input.Contains("plex", StringComparison.OrdinalIgnoreCase))
        {
            yield return input.Replace("plex", "-plex", StringComparison.OrdinalIgnoreCase);
            yield return input.Replace("plex", "", StringComparison.OrdinalIgnoreCase);
        }

        if (input.Contains("TMT18", StringComparison.OrdinalIgnoreCase))
        {
            yield return input.Replace("TMT18", "TMTpro", StringComparison.OrdinalIgnoreCase);
        }
        else if (input.Contains("TMTpro", StringComparison.OrdinalIgnoreCase))
        {
            yield return input.Replace("TMTpro", "TMT18", StringComparison.OrdinalIgnoreCase);
        }

        if (input.Contains("TMT", StringComparison.OrdinalIgnoreCase) && input.Contains("pro", StringComparison.OrdinalIgnoreCase))
        {
            yield return input.Replace("pro", "", StringComparison.OrdinalIgnoreCase);
        }
    }

    #endregion

    #region Utilities

    protected virtual bool MatchesIdentifier(Modification modification, string identifier)
    {
        if (string.IsNullOrWhiteSpace(identifier))
        {
            return false;
        }

        // Direct match on IdWithMotif
        if (string.Equals(modification.IdWithMotif, identifier, StringComparison.OrdinalIgnoreCase))
        {
            return true;
        }

        // Match on OriginalId
        if (string.Equals(modification.OriginalId, identifier, StringComparison.OrdinalIgnoreCase))
        {
            return true;
        }

        // Match on Type:IdWithMotif pattern
        if (modification is { ModificationType: not null, IdWithMotif: not null })
        {
            var fullId = $"{modification.ModificationType}:{modification.IdWithMotif}";
            if (string.Equals(fullId, identifier, StringComparison.OrdinalIgnoreCase))
            {
                return true;
            }
        }

        // If identifier contains a colon, try extracting just the suffix (after colon)
        var colonIndex = identifier.IndexOf(':');
        if (colonIndex > 0 && colonIndex < identifier.Length - 1)
        {
            var suffix = identifier[(colonIndex + 1)..];
            if (string.Equals(modification.IdWithMotif, suffix, StringComparison.OrdinalIgnoreCase) ||
                string.Equals(modification.OriginalId, suffix, StringComparison.OrdinalIgnoreCase))
            {
                return true;
            }
        }

        return false;
    }

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

    /// <summary>
    /// Checks if the modification is already resolved to this lookup's database.
    /// Used to avoid redundant lookups when the mod is already in the correct format.
    /// </summary>
    protected virtual bool IsResolvedToThisDatabase(CanonicalModification mod)
    {
        if (!mod.IsResolved || mod.MzLibModification == null)
        {
            return false;
        }

        // Check if the modification type matches this lookup's database name
        return string.Equals(mod.MzLibModification.ModificationType, Name, StringComparison.OrdinalIgnoreCase);
    }

    #endregion
}
