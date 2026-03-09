using Chemistry;
using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Extension methods for <see cref="CanonicalSequenceBuilder"/> to provide
/// convenient overloads for adding modifications without requiring originalRepresentation.
/// </summary>
public static class CanonicalSequenceBuilderExtensions
{
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
