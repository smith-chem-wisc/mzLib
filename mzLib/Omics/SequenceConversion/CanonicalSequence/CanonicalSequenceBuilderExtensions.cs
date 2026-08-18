using Chemistry;
using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Extension methods for <see cref="CanonicalSequenceBuilder"/> to provide
/// convenient overloads for adding modifications without requiring originalRepresentation,
/// and for converting from mzLib biopolymer objects.
/// </summary>
public static class CanonicalSequenceBuilderExtensions
{
    /// <summary>
    /// Creates a <see cref="CanonicalSequenceBuilder"/> from an <see cref="IBioPolymerWithSetMods"/>.
    /// Converts the biopolymer's modifications to canonical format.
    /// 
    /// IMPORTANT: AllModsOneIsNterminus indexing convention:
    /// - Index 1: N-terminal modification
    /// - Index r+2: Residue modification at zero-based position r
    /// - Index Length+2: C-terminal modification
    /// </summary>
    /// <param name="bioPolymer">The biopolymer with modifications to convert.</param>
    /// <returns>A builder populated with the biopolymer's sequence and modifications.</returns>
    public static CanonicalSequenceBuilder ToCanonicalSequenceBuilder(this IBioPolymerWithSetMods bioPolymer)
    {
        var builder = new CanonicalSequenceBuilder(bioPolymer.BaseSequence)
            .WithSourceFormat("mzLib");

        // Add N-terminal modification if present
        // AllModsOneIsNterminus uses index 1 for N-terminus
        if (bioPolymer.AllModsOneIsNterminus.TryGetValue(1, out Modification? nTermMod))
        {
            builder.AddNTerminalModification(
                originalRepresentation: nTermMod.IdWithMotif,
                mass: nTermMod.MonoisotopicMass,
                formula: nTermMod.ChemicalFormula,
                mzLibId: nTermMod.IdWithMotif,
                mzLibModification: nTermMod);
        }

        // Add residue modifications
        // AllModsOneIsNterminus uses index r+2 for residue at zero-based index r
        for (int r = 0; r < bioPolymer.Length; r++)
        {
            if (bioPolymer.AllModsOneIsNterminus.TryGetValue(r + 2, out Modification? residueMod))
            {
                builder.AddResidueModification(
                    residueIndex: r,
                    originalRepresentation: residueMod.IdWithMotif,
                    mass: residueMod.MonoisotopicMass,
                    formula: residueMod.ChemicalFormula,
                    mzLibId: residueMod.IdWithMotif,
                    mzLibModification: residueMod);
            }
        }

        // Add C-terminal modification if present
        // AllModsOneIsNterminus uses index Length+2 for C-terminus
        if (bioPolymer.AllModsOneIsNterminus.TryGetValue(bioPolymer.Length + 2, out Modification? cTermMod))
        {
            builder.AddCTerminalModification(
                originalRepresentation: cTermMod.IdWithMotif,
                mass: cTermMod.MonoisotopicMass,
                formula: cTermMod.ChemicalFormula,
                mzLibId: cTermMod.IdWithMotif,
                mzLibModification: cTermMod);
        }

        return builder;
    }
    /// <summary>
    /// Adds a residue modification at the specified zero-based index without requiring originalRepresentation.
    /// The mzLibId or mass will be used as the original representation.
    /// </summary>
    public static CanonicalSequenceBuilder AddResidueModification(
        this CanonicalSequenceBuilder builder,
        int residueIndex,
        double? mass = null,
        ChemicalFormula? formula = null,
        int? unimodId = null,
        string? mzLibId = null,
        Modification? mzLibModification = null)
    {
        // Use mzLibId as originalRepresentation if available, otherwise use mass string
        string originalRep = mzLibId ?? (mass.HasValue ? $"+{mass.Value}" : "Unknown");
        
        return builder.AddResidueModification(
            residueIndex, originalRep, mass, formula, unimodId, mzLibId, mzLibModification);
    }

    /// <summary>
    /// Adds an N-terminal modification without requiring originalRepresentation.
    /// The mzLibId or mass will be used as the original representation.
    /// </summary>
    public static CanonicalSequenceBuilder AddNTerminalModification(
        this CanonicalSequenceBuilder builder,
        double? mass = null,
        ChemicalFormula? formula = null,
        int? unimodId = null,
        string? mzLibId = null,
        Modification? mzLibModification = null)
    {
        // Use mzLibId as originalRepresentation if available, otherwise use mass string
        string originalRep = mzLibId ?? (mass.HasValue ? $"+{mass.Value}" : "Unknown");
        
        return builder.AddNTerminalModification(
            originalRep, mass, formula, unimodId, mzLibId, mzLibModification);
    }

    /// <summary>
    /// Adds a C-terminal modification without requiring originalRepresentation.
    /// The mzLibId or mass will be used as the original representation.
    /// </summary>
    public static CanonicalSequenceBuilder AddCTerminalModification(
        this CanonicalSequenceBuilder builder,
        double? mass = null,
        ChemicalFormula? formula = null,
        int? unimodId = null,
        string? mzLibId = null,
        Modification? mzLibModification = null)
    {
        // Use mzLibId as originalRepresentation if available, otherwise use mass string
        string originalRep = mzLibId ?? (mass.HasValue ? $"+{mass.Value}" : "Unknown");
        
        return builder.AddCTerminalModification(
            originalRep, mass, formula, unimodId, mzLibId, mzLibModification);
    }
}
