namespace Omics.BioPolymerGroup;

/// <summary>
/// Identifies the type of biopolymer in a <see cref="BioPolymerGroup"/>,
/// which primarily determines the occupancy calculation strategy.
/// </summary>
public enum BioPolymerGroupType
{
    /// <summary>
    /// Protein group — occupancy is calculated at protein-level coordinates
    /// using <see cref="ModificationOccupancyCalculator.CalculateProteinLevelOccupancy"/>.
    /// </summary>
    Protein,

    /// <summary>
    /// Peptide group — occupancy is calculated in peptide-local coordinates
    /// using <see cref="ModificationOccupancyCalculator.CalculatePeptideLevelOccupancy"/>.
    /// </summary>
    Peptide,

    /// <summary>
    /// Oligonucleotide group — occupancy is calculated in oligo-local coordinates
    /// using <see cref="ModificationOccupancyCalculator.CalculatePeptideLevelOccupancy"/>.
    /// </summary>
    Oligo
}
