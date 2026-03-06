using Chemistry;
using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Specifies the position type of a modification within a sequence.
/// </summary>
public enum ModificationPositionType
{
    /// <summary>
    /// Modification is on the N-terminus (or 5' end for nucleic acids).
    /// </summary>
    NTerminus,

    /// <summary>
    /// Modification is on the C-terminus (or 3' end for nucleic acids).
    /// </summary>
    CTerminus,

    /// <summary>
    /// Modification is on a specific residue within the sequence.
    /// </summary>
    Residue
}

/// <summary>
/// Immutable representation of a modification in the canonical sequence format.
/// This serves as the universal intermediate representation during sequence conversions.
/// Fields are populated lazily from the source format - only what the source provides
/// is guaranteed to be present. Enrichment can add additional fields during serialization.
/// </summary>
/// <param name="PositionType">Whether this mod is on N-terminus, C-terminus, or a residue.</param>
/// <param name="ResidueIndex">Zero-based index of the modified residue in the base sequence. 
/// Null for terminal modifications.</param>
/// <param name="TargetResidue">The amino acid or nucleotide this modification is attached to.
/// Null for terminal modifications that don't target a specific residue.</param>
/// <param name="OriginalRepresentation">The raw string representation from the source format.
/// Always populated - this is what was parsed from the input.</param>
/// <param name="MonoisotopicMass">The monoisotopic mass shift of this modification.
/// May be null if the source format doesn't provide mass information.</param>
/// <param name="ChemicalFormula">The chemical formula of this modification.
/// May be null if the source format doesn't provide formula information.</param>
/// <param name="UnimodId">The UNIMOD accession number (e.g., 35 for Oxidation).
/// May be null if the source format doesn't use UNIMOD identifiers.</param>
/// <param name="MzLibId">The mzLib modification identifier (e.g., "Oxidation on M").
/// May be null if the source format doesn't use mzLib identifiers.</param>
/// <param name="MzLibModification">Reference to the resolved mzLib Modification object.
/// May be null if not yet resolved or if no matching modification exists.</param>
public readonly record struct CanonicalModification(
    ModificationPositionType PositionType,
    int? ResidueIndex,
    char? TargetResidue,
    string OriginalRepresentation,
    double? MonoisotopicMass = null,
    ChemicalFormula? ChemicalFormula = null,
    int? UnimodId = null,
    string? MzLibId = null,
    Modification? MzLibModification = null)
{
    /// <summary>
    /// Creates a residue modification at the specified zero-based index.
    /// </summary>
    public static CanonicalModification AtResidue(
        int residueIndex,
        char targetResidue,
        string originalRepresentation,
        double? mass = null,
        ChemicalFormula? formula = null,
        int? unimodId = null,
        string? mzLibId = null,
        Modification? mzLibModification = null)
    {
        return new CanonicalModification(
            ModificationPositionType.Residue,
            residueIndex,
            targetResidue,
            originalRepresentation,
            mass,
            formula,
            unimodId,
            mzLibId,
            mzLibModification);
    }

    /// <summary>
    /// Creates an N-terminal (or 5') modification.
    /// </summary>
    public static CanonicalModification AtNTerminus(
        string originalRepresentation,
        char? targetResidue = null,
        double? mass = null,
        ChemicalFormula? formula = null,
        int? unimodId = null,
        string? mzLibId = null,
        Modification? mzLibModification = null)
    {
        return new CanonicalModification(
            ModificationPositionType.NTerminus,
            null,
            targetResidue,
            originalRepresentation,
            mass,
            formula,
            unimodId,
            mzLibId,
            mzLibModification);
    }

    /// <summary>
    /// Creates a C-terminal (or 3') modification.
    /// </summary>
    public static CanonicalModification AtCTerminus(
        string originalRepresentation,
        char? targetResidue = null,
        double? mass = null,
        ChemicalFormula? formula = null,
        int? unimodId = null,
        string? mzLibId = null,
        Modification? mzLibModification = null)
    {
        return new CanonicalModification(
            ModificationPositionType.CTerminus,
            null,
            targetResidue,
            originalRepresentation,
            mass,
            formula,
            unimodId,
            mzLibId,
            mzLibModification);
    }

    /// <summary>
    /// Returns true if this modification has mass information available.
    /// </summary>
    public bool HasMass => MonoisotopicMass.HasValue || ChemicalFormula != null;

    /// <summary>
    /// Gets the effective mass of this modification, calculating from formula if necessary.
    /// Returns null if no mass information is available.
    /// </summary>
    public double? EffectiveMass => MonoisotopicMass ?? ChemicalFormula?.MonoisotopicMass ?? MzLibModification?.MonoisotopicMass;

    /// <summary>
    /// Returns true if this modification has been resolved to an mzLib Modification.
    /// </summary>
    public bool IsResolved => MzLibModification != null;

    /// <summary>
    /// Creates a new CanonicalModification with the specified mzLib Modification resolved.
    /// This is used during enrichment to populate additional fields.
    /// </summary>
    public CanonicalModification WithResolvedModification(Modification modification)
    {
        return this with
        {
            MzLibModification = modification,
            MzLibId = modification.IdWithMotif,
            MonoisotopicMass = MonoisotopicMass ?? modification.MonoisotopicMass,
            ChemicalFormula = ChemicalFormula ?? modification.ChemicalFormula
        };
    }

    /// <summary>
    /// Creates a new CanonicalModification with the specified mass.
    /// </summary>
    public CanonicalModification WithMass(double mass)
    {
        return this with { MonoisotopicMass = mass };
    }

    public override string ToString()
    {
        var position = PositionType switch
        {
            ModificationPositionType.NTerminus => "N-term",
            ModificationPositionType.CTerminus => "C-term",
            ModificationPositionType.Residue => $"@{ResidueIndex}({TargetResidue})",
            _ => "Unknown"
        };

        var identifier = MzLibId ?? (UnimodId.HasValue ? $"UNIMOD:{UnimodId}" : OriginalRepresentation);
        return $"{identifier} {position}";
    }
}
