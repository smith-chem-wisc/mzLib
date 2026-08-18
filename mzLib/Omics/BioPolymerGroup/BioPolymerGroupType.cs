namespace Omics.BioPolymerGroup;

/// <summary>
/// Identifies the type of biopolymer in a <see cref="BioPolymerGroup"/>,
/// which primarily determines the occupancy calculation position-mapping strategy.
/// </summary>
public enum BioPolymerGroupType
{
    /// <summary>
    /// Parent (Protein/NucleicAcid) group — occupancy is calculated at protein-level positions
    /// using <see cref="ModificationOccupancyCalculator.CalculateParentLevelOccupancy"/>.
    /// </summary>
    Parent,

    /// <summary>
    /// DigestionProduct (Peptide/Oligo) group — occupancy is calculated in product-local positions
    /// using <see cref="ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy"/>.
    /// </summary>
    DigestionProduct
}
