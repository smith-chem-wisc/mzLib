namespace Omics.Modifications.IO.Modomics;

/// <summary>
/// Why a MODOMICS entry could not be represented as an mzLib residue modification.
/// </summary>
public enum ModomicsRepresentationFailureReason
{
    /// <summary>The source entry has no chemical formula.</summary>
    EmptyFormula,

    /// <summary>The entry is a stand-in for an unknown modification.</summary>
    UnknownEntry,

    /// <summary>The source formula could not be parsed.</summary>
    InvalidFormula,

    /// <summary>The entry targets a generic moiety (e.g. "X") rather than a specific base.</summary>
    UnsupportedReferenceMoiety,

    /// <summary>The CSV moiety type is not one of nucleoside/nucleotide/base.</summary>
    UnsupportedMoietyType,

    /// <summary>A modified base (e.g. preQ0) requires a new residue definition, not a mass shift.</summary>
    BaseMoietyRequiresNewResidue,

    /// <summary>
    /// Real-base nucleotide formulas are published as inconsistent phosphate anions (e.g. pm1G is
    /// "C11H14N5O8P", two protons shy of the free acid), so an exact neutral shift cannot be derived;
    /// the neutral nucleoside entries carry the same chemistry.
    /// </summary>
    NucleotideProtonationAmbiguous,

    /// <summary>
    /// The entry describes a terminus state or terminal cofactor (e.g. "5' diphosphate end", NAD caps)
    /// that belongs to the five/three prime terminus model rather than a residue mass shift.
    /// </summary>
    TerminalStructureRequiresTerminusModel,

    /// <summary>
    /// The conversion yields an empty formula: the entry imparts no mass shift (canonical nucleosides,
    /// isomerizations such as pseudouridine) and cannot be represented as a mass-difference modification.
    /// </summary>
    NoMassShift,

    /// <summary>The reference moiety does not map to a valid modification motif.</summary>
    InvalidTargetMotif
}
