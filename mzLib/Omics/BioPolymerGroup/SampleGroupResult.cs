using MassSpectrometry;

namespace Omics.BioPolymerGroup;

/// <summary>
/// Bundles quantification and modification occupancy data for a single sample group
/// (Condition × BiologicalReplicate). Each group contributes 2 columns (SpectralCount + CountOccupancy)
/// or 4 columns (+Intensity + IntensityOccupancy) when intensity data is available.
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
    /// The sample files (or channels) that belong to this result group.
    /// For label-free data, contains one or more <see cref="SpectraFileInfo"/> entries (one per fraction).
    /// For isobaric data, contains a single <see cref="IsobaricQuantSampleInfo"/> entry per channel.
    /// </summary>
    public Dictionary<string, ISampleInfo> FilesInGroup { get; init; } = new();

    /// <summary>
    /// Per-file intensity values for this result group, keyed by sample info.
    /// Populated from <see cref="BioPolymerGroup.IntensitiesBySample"/> filtered to the files in this group.
    /// </summary>
    public Dictionary<string, double> IntensitiesBySample { get; init; } = new();

    /// <summary>
    /// Summed intensity across all files in this sample group.
    /// Computed from <see cref="IntensitiesBySample"/>. Zero when no intensity data is available.
    /// </summary>
    public double Intensity => IntensitiesBySample.Values.Sum();

    /// <summary>
    /// True when intensity data was available for this sample group (i.e., <see cref="IntensitiesBySample"/> is non-empty).
    /// Controls whether intensity and intensity-occupancy columns are output.
    /// </summary>
    public bool HasIntensityData => IntensitiesBySample.Count > 0;

    #endregion

    #region Occupancy

    /// <summary>
    /// Protein-level modification occupancy keyed by biopolymer accession, then by one-based protein position.
    /// Populated by <see cref="ModificationOccupancyCalculator.CalculateParentLevelOccupancy"/>.
    /// </summary>
    public Dictionary<string, Dictionary<int, List<SiteSpecificModificationOccupancy>>> ParentOccupancy { get; } = new();

    /// <summary>
    /// Peptide-level modification occupancy keyed by base sequence, then by peptide-local position
    /// (AllModsOneIsNterminus convention: 1 = N-terminus, 2 = first residue, etc.).
    /// Populated by <see cref="ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy"/>.
    /// </summary>
    public Dictionary<string, Dictionary<int, List<SiteSpecificModificationOccupancy>>> DigestionProductOccupancy { get; } = new();

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
        var occupancy = proteinLevel ? ParentOccupancy : DigestionProductOccupancy;
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
        var occupancy = proteinLevel ? ParentOccupancy : DigestionProductOccupancy;
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