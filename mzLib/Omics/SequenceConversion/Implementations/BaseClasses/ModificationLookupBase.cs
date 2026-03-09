using Chemistry;
using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Base class for modification lookups providing shared resolution logic.
/// Subclasses implement database-specific lookup strategies.
/// </summary>
public abstract class ModificationLookupBase : IModificationLookup
{
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
    /// </summary>
    /// <param name="name">The modification name or identifier.</param>
    /// <param name="targetResidue">Optional target residue for disambiguation.</param>
    /// <returns>The resolved Modification, or null if not found.</returns>
    protected abstract Modification? TryResolveByName(string name, char? targetResidue);

    /// <summary>
    /// Attempts to find a modification matching the given chemical formula.
    /// Default implementation returns null; subclasses can override.
    /// </summary>
    /// <param name="formula">The chemical formula to match.</param>
    /// <param name="targetResidue">Optional target residue for disambiguation.</param>
    /// <returns>The resolved Modification, or null if not found.</returns>
    protected virtual Modification? TryResolveByFormula(ChemicalFormula formula, char? targetResidue) => null;

    /// <summary>
    /// Attempts to find a modification matching the given mass within tolerance.
    /// Default implementation returns null; subclasses can override.
    /// </summary>
    /// <param name="mass">The monoisotopic mass to match.</param>
    /// <param name="targetResidue">Optional target residue for disambiguation.</param>
    /// <returns>The resolved Modification, or null if not found.</returns>
    protected virtual Modification? TryResolveByMass(double mass, char? targetResidue) => null;

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
    /// Helper method to select a modification from candidates with residue preference.
    /// Prefers modifications that target the specified residue if provided.
    /// </summary>
    /// <param name="candidates">The candidate modifications to select from.</param>
    /// <param name="targetResidue">Optional preferred target residue.</param>
    /// <returns>The best matching modification, or null if no candidates.</returns>
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
}
