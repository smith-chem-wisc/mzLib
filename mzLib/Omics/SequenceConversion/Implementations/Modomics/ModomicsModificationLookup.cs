using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Resolves MODOMICS one-letter codes and identifiers to MODOMICS RNA
/// modifications, keeping them separate from chemistry-equivalent curated mods.
/// </summary>
public sealed class ModomicsModificationLookup : ModificationLookupBase
{
    public static ModomicsModificationLookup Instance { get; } = new();

    private ModomicsModificationLookup()
        : base(Mods.ModomicsRnaModifications, null)
    {
    }

    public override string Name => "Modomics";

    public override CanonicalModification? TryResolve(
        string originalRepresentation,
        char? targetResidue = null,
        Chemistry.ChemicalFormula? chemicalFormula = null,
        ModificationPositionType? positionType = null)
    {
        if (originalRepresentation?.Length == 1 &&
            TryResolveCode(originalRepresentation[0], targetResidue, out var modification))
        {
            return new CanonicalModification(
                positionType ?? ModificationPositionType.Residue,
                null,
                targetResidue ?? GetTargetResidue(modification),
                originalRepresentation,
                modification.MonoisotopicMass,
                modification.ChemicalFormula,
                MzLibId: modification.IdWithMotif,
                MzLibModification: modification);
        }

        return base.TryResolve(originalRepresentation, targetResidue, chemicalFormula, positionType);
    }

    public bool TryResolveCode(char code, char? targetResidue, out Modification? modification)
    {
        modification = null;
        if (!Mods.ModomicsLoadReport.ModificationsByAbbreviation.TryGetValue(code.ToString(), out var candidates))
        {
            return false;
        }

        var matchingCandidates = candidates
            .Where(candidate => targetResidue is null || string.Equals(
                GetTargetResidue(candidate).ToString(), targetResidue.Value.ToString(), StringComparison.OrdinalIgnoreCase))
            .ToList();

        if (matchingCandidates.Count != 1)
        {
            return false;
        }

        modification = matchingCandidates[0];
        return true;
    }

    private static char GetTargetResidue(Modification modification) =>
        modification.Target.Motif?.FirstOrDefault() ?? '\0';
}
