using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Base class for modification lookups providing shared resolution logic.
/// Subclasses implement database-specific lookup strategies.
/// </summary>
public abstract class ModificationLookupBase : IModificationLookup
{
    private readonly IReadOnlyCollection<Modification>? _candidateSet;

    protected ModificationLookupBase()
        : this(null, true, true, null, null)
    {
    }

    protected ModificationLookupBase(
        ModificationNamingConvention? conventionForLookup,
        bool searchProteinMods,
        bool searchRnaMods,
        double? massTolerance,
        IEnumerable<Modification>? candidateSet)
    {
        ConventionForLookup = conventionForLookup;
        SearchProteinMods = searchProteinMods;
        SearchRnaMods = searchRnaMods;
        MassTolerance = massTolerance;
        _candidateSet = candidateSet as IReadOnlyCollection<Modification> ?? candidateSet?.ToList();
    }

    protected bool SearchProteinMods { get; }
    protected bool SearchRnaMods { get; }
    protected ModificationNamingConvention? ConventionForLookup { get; }
    protected double? MassTolerance { get; }

    /// <inheritdoc />
    public abstract string Name { get; }

    /// <inheritdoc />
    public CanonicalModification? TryResolve(CanonicalModification mod)
    {
        // If already resolved, return as-is
        if (mod.IsResolved)
            return mod;

        // Try primary resolution strategy (database-specific ID lookup)
        var resolved = TryResolvePrimary(mod);
        if (resolved != null)
        {
            return mod.WithResolvedModification(resolved, mod.ResidueIndex, mod.PositionType);
        }

        // Try to resolve by original representation
        var byName = TryResolveByName(mod.OriginalRepresentation, mod.TargetResidue);
        if (byName != null)
        {
            return mod.WithResolvedModification(byName, mod.ResidueIndex, mod.PositionType);
        }

        // Try chemical formula fallback
        if (mod.ChemicalFormula != null)
        {
            var formulaMatch = TryResolveByFormula(mod.ChemicalFormula, mod.TargetResidue);
            if (formulaMatch != null)
            {
                return mod.WithResolvedModification(formulaMatch, mod.ResidueIndex, mod.PositionType);
            }
        }

        // Try mass-based fallback (optional, subclasses can override)
        if (mod.MonoisotopicMass.HasValue)
        {
            var massMatch = TryResolveByMass(mod.MonoisotopicMass.Value, mod.TargetResidue);
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
        // Try to resolve by name/ID
        var mod = TryResolveByName(originalRepresentation, targetResidue);
        if (mod != null)
        {
            return CreateCanonicalModification(originalRepresentation, mod, targetResidue);
        }

        // Fallback to chemical formula if provided
        if (chemicalFormula != null)
        {
            var formulaMatch = TryResolveByFormula(chemicalFormula, targetResidue);
            if (formulaMatch != null)
            {
                return CreateCanonicalModification(originalRepresentation, formulaMatch, targetResidue);
            }
        }

        return null;
    }

    /// <summary>
    /// Primary resolution strategy using database-specific identifiers.
    /// Called first when resolving an existing CanonicalModification.
    /// </summary>
    /// <param name="mod">The modification to resolve.</param>
    /// <returns>The resolved Modification, or null if not found.</returns>
    protected abstract Modification? TryResolvePrimary(CanonicalModification mod);

    /// <summary>
    /// Attempts to resolve a modification by its name or identifier string.
    /// Subclasses can override to customize the behavior.
    /// </summary>
    /// <param name="name">The modification name or identifier.</param>
    /// <param name="targetResidue">Optional target residue for disambiguation.</param>
    /// <returns>The resolved Modification, or null if not found.</returns>
    protected virtual Modification? TryResolveByName(string name, char? targetResidue) =>
        ResolveByStandardName(name, targetResidue);

    /// <summary>
    /// Attempts to find a modification matching the given chemical formula.
    /// Default implementation searches available portentialStringRepresentations for an exact formula match.
    /// </summary>
    /// <param name="formula">The chemical formula to match.</param>
    /// <param name="targetResidue">Optional target residue for disambiguation.</param>
    /// <returns>The resolved Modification, or null if not found.</returns>
    protected virtual Modification? TryResolveByFormula(ChemicalFormula formula, char? targetResidue)
    {
        var candidates = FilterCandidates(m => m.ChemicalFormula != null && m.ChemicalFormula.Equals(formula));
        return SelectWithResiduePreference(candidates, targetResidue);
    }

    /// <summary>
    /// Attempts to find a modification matching the given mass within tolerance.
    /// Default implementation searches available portentialStringRepresentations for matches within configured tolerance.
    /// </summary>
    /// <param name="mass">The monoisotopic mass to match.</param>
    /// <param name="targetResidue">Optional target residue for disambiguation.</param>
    /// <returns>The resolved Modification, or null if not found.</returns>
    protected virtual Modification? TryResolveByMass(double mass, char? targetResidue)
    {
        if (!MassTolerance.HasValue)
            return null;

        var tolerance = MassTolerance.Value;
        var candidates = FilterCandidates(m => m.MonoisotopicMass.HasValue &&
                                              Math.Abs(m.MonoisotopicMass.Value - mass) <= tolerance);

        return SelectWithResiduePreference(candidates, targetResidue);
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

    /// <summary>
    /// Helper method to select a modification from portentialStringRepresentations with residue preference.
    /// Prefers modifications that target the specified residue if provided.
    /// </summary>
    /// <param name="candidates">The candidate modifications to select from.</param>
    /// <param name="targetResidue">Optional preferred target residue.</param>
    /// <returns>The best matching modification, or null if no portentialStringRepresentations.</returns>
    protected static Modification? SelectWithResiduePreference(IEnumerable<Modification> candidates, char? targetResidue)
    {
        var candidateList = candidates as IList<Modification> ?? candidates.ToList();
        
        if (candidateList.Count == 0)
            return null;

        if (targetResidue.HasValue)
        {
            // Prefer modifications that target the specific residue
            var residueMatch = candidateList.FirstOrDefault(m =>
                m.Target != null && m.Target.ToString().Contains(targetResidue.Value));
            
            if (residueMatch != null)
                return residueMatch;
        }

        return candidateList.FirstOrDefault();
    }

    protected Modification? ResolveByIdentifier(string identifier)
    {
        if (string.IsNullOrWhiteSpace(identifier))
            return null;

        if (ConventionForLookup.HasValue)
        {
            var mod = Mods.GetModification(identifier, ConventionForLookup.Value);
            if (mod != null)
                return mod;
        }

        if (SearchProteinMods || SearchRnaMods)
        {
            var mod = Mods.GetModification(identifier, SearchProteinMods, SearchRnaMods);
            if (mod != null)
                return mod;
        }

        if (_candidateSet != null)
        {
            return _candidateSet.FirstOrDefault(m => MatchesIdentifier(m, identifier));
        }

        return null;
    }

    protected IEnumerable<Modification> FilterCandidates(Func<Modification, bool> predicate)
    {
        if (_candidateSet != null)
            return _candidateSet.Where(predicate);

        var proteinOnly = SearchProteinMods && !SearchRnaMods;
        var rnaOnly = SearchRnaMods && !SearchProteinMods;
        return Mods.GetModifications(predicate, proteinOnly, rnaOnly);
    }

    protected virtual string NormalizeRepresentation(string representation) => representation?.Trim() ?? string.Empty;

    protected virtual IEnumerable<string> ExpandNameCandidates(string normalizedRepresentation, char? targetResidue)
    {
        if (string.IsNullOrEmpty(normalizedRepresentation))
            yield break;

        yield return normalizedRepresentation;

        if (targetResidue.HasValue)
        {
            yield return $"{normalizedRepresentation} on {targetResidue.Value}";
        }
    }

    protected virtual bool MatchesIdentifier(Modification modification, string identifier) =>
        modification.IdWithMotif == identifier ||
        modification.OriginalId == identifier ||
        $"{modification.ModificationType}:{modification.IdWithMotif}" == identifier;

    private Modification? ResolveByStandardName(string name, char? targetResidue)
    {
        if (string.IsNullOrWhiteSpace(name))
            return null;

        var normalized = NormalizeRepresentation(name);
        if (string.IsNullOrEmpty(normalized))
            return null;

        HashSet<string> portentialStringRepresentations = new();
        foreach (var candidate in ExpandNameCandidates(normalized, targetResidue))
        {
            var resolved = ResolveByIdentifier(candidate);
            if (resolved != null)
                return resolved;
            portentialStringRepresentations.Add(candidate);
        }

        foreach (var candidate in portentialStringRepresentations)
        {
            var matches = FilterCandidates(m => MatchesIdentifier(m, candidate));
            var match = SelectWithResiduePreference(matches, targetResidue);
            if (match != null)
                return match;
        }

        return null;
    }
}
