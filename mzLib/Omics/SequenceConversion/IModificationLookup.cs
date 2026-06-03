using Chemistry;

namespace Omics.SequenceConversion;

/// <summary>
/// Resolves and enriches modifications by looking up additional information
/// from a modification database or dictionary.
/// 
/// Implementations may resolve by:
/// - UNIMOD ID (e.g., "UNIMOD:35" → Oxidation)
/// - mzLib ID (e.g., "Oxidation on M")
/// - Mass (fuzzy matching within tolerance)
/// - Chemical formula
/// - Name patterns
/// </summary>
public interface IModificationLookup
{
    /// <summary>
    /// Gets a descriptive name for this lookup (e.g., "UNIMOD", "mzLib", "MassBased").
    /// </summary>
    string Name { get; }

    /// <summary>
    /// Attempts to resolve a modification by looking up additional information.
    /// Returns an enriched modification if found, or null if not found.
    /// 
    /// The returned modification will have additional fields populated
    /// (e.g., MzLibModification, ChemicalFormula, MonoisotopicMass) based on
    /// what the lookup source provides.
    /// </summary>
    /// <param name="mod">The modification to resolve.</param>
    /// <returns>An enriched modification if found; otherwise, null.</returns>
    CanonicalModification? TryResolve(CanonicalModification mod);

    /// <summary>
    /// Attempts to resolve a modification by its original string representation.
    /// This is useful when parsing and you only have the raw string.
    /// </summary>
    /// <param name="originalRepresentation">The raw modification string from the source format.</param>
    /// <param name="targetResidue">The amino acid or nucleotide this mod is attached to, if known.</param>
    /// <param name="chemicalFormula">Optional chemical formula to use as a fallback for matching
    /// when the original representation cannot be resolved directly.</param>
    /// <param name="positionType">N, C, or internal</param>
    /// <returns>An enriched modification if found; otherwise, null.</returns>
    CanonicalModification? TryResolve(string originalRepresentation, char? targetResidue = null, ChemicalFormula? chemicalFormula = null, ModificationPositionType? positionType = null);
}
