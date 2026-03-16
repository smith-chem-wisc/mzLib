namespace Omics.BioPolymerGroup;

/// <summary>
/// Bundles quantification and modification occupancy data for a single sample group
/// (Condition × BiologicalReplicate). Each group contributes 2 columns (SpectralCount + CountOccupancy)
/// or 4 columns (+Intensity + IntensityOccupancy) when FlashLFQ intensity data is available.
/// </summary>
public sealed class SampleGroupResult
{
    #region Identity

    /// <summary>
    /// Experimental condition (e.g., "Control", "Treatment"). May be empty for simple designs.
    /// </summary>
    public string Condition { get; }

    /// <summary>
    /// Biological replicate index within the condition.
    /// </summary>
    public int BiologicalReplicate { get; }

    /// <summary>
    /// Display label for column headers (e.g., "Control_1" or a filename).
    /// Set by the caller based on experimental design context.
    /// </summary>
    public string Label { get; init; } = string.Empty;

    #endregion

    #region Quantification

    /// <summary>
    /// Number of PSMs (spectral matches) in this sample group for this biopolymer group.
    /// </summary>
    public int SpectralCount { get; set; }

    /// <summary>
    /// Summed intensity from FlashLFQ for this sample group. Only meaningful when <see cref="HasIntensityData"/> is true.
    /// </summary>
    public double Intensity { get; set; }

    /// <summary>
    /// True when FlashLFQ intensity data was available for this sample group.
    /// Controls whether intensity and intensity-occupancy columns are output.
    /// </summary>
    public bool HasIntensityData { get; init; }

    #endregion

    #region Occupancy

    /// <summary>
    /// Protein-level modification occupancy keyed by biopolymer accession, then by one-based protein position.
    /// Populated by <see cref="ModificationOccupancyCalculator.CalculateProteinLevelOccupancy"/>.
    /// </summary>
    public Dictionary<string, Dictionary<int, List<SiteSpecificModificationOccupancy>>> ProteinOccupancy { get; } = new();

    /// <summary>
    /// Peptide-level modification occupancy keyed by base sequence, then by peptide-local position
    /// (AllModsOneIsNterminus convention: 1 = N-terminus, 2 = first residue, etc.).
    /// Populated by <see cref="ModificationOccupancyCalculator.CalculatePeptideLevelOccupancy"/>.
    /// </summary>
    public Dictionary<string, Dictionary<int, List<SiteSpecificModificationOccupancy>>> PeptideOccupancy { get; } = new();

    #endregion

    public SampleGroupResult(string condition, int biologicalReplicate)
    {
        Condition = condition;
        BiologicalReplicate = biologicalReplicate;
    }

    #region Formatting

    /// <summary>
    /// Formats count-based occupancy for a TSV cell.
    /// Output: semicolon-separated mod entries within each entity, pipe-separated between entities.
    /// </summary>
    /// <param name="orderedKeys">Ordered accessions (protein-level) or base sequences (peptide-level).</param>
    /// <param name="proteinLevel">True for protein-level occupancy; false for peptide-level.</param>
    public string FormatCountOccupancy(IEnumerable<string> orderedKeys, bool proteinLevel = true)
    {
        var occupancy = proteinLevel ? ProteinOccupancy : PeptideOccupancy;
        return FormatOccupancy(occupancy, orderedKeys, o => o.ToSpectralCountModInfoString());
    }

    /// <summary>
    /// Formats intensity-based stoichiometry for a TSV cell.
    /// Only meaningful when <see cref="HasIntensityData"/> is true.
    /// </summary>
    /// <param name="orderedKeys">Ordered accessions (protein-level) or base sequences (peptide-level).</param>
    /// <param name="proteinLevel">True for protein-level occupancy; false for peptide-level.</param>
    public string FormatIntensityOccupancy(IEnumerable<string> orderedKeys, bool proteinLevel = true)
    {
        var occupancy = proteinLevel ? ProteinOccupancy : PeptideOccupancy;
        return FormatOccupancy(occupancy, orderedKeys, o => o.ToIntensityModInfoString());
    }

    /// <summary>
    /// Core formatting helper. Iterates ordered keys, formats each entity's modifications,
    /// and joins with the standard separators (; within entity, | between entities).
    /// </summary>
    private static string FormatOccupancy(
        Dictionary<string, Dictionary<int, List<SiteSpecificModificationOccupancy>>> occupancy,
        IEnumerable<string> orderedKeys,
        Func<SiteSpecificModificationOccupancy, string> formatter)
    {
        var parts = new List<string>();

        foreach (var key in orderedKeys)
        {
            if (!occupancy.TryGetValue(key, out var positions))
                continue;

            string entityString = string.Join(";",
                positions.OrderBy(kvp => kvp.Key)
                    .SelectMany(kvp => kvp.Value)
                    .Select(formatter));

            if (!string.IsNullOrEmpty(entityString))
                parts.Add(entityString);
        }

        return string.Join("|", parts);
    }

    #endregion
}