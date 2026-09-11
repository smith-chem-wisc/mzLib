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
    /// <remarks>
    /// Not unique. Whenever the label is derived from a file name, two files with the same name in
    /// different directories produce the same label. Use <see cref="Identity"/> to tell sample
    /// groups apart; a label is for display only.
    /// </remarks>
    public string Label { get; init; } = string.Empty;

    /// <summary>
    /// Stable identifier for this sample group, unique within a dataset. What identifies a sample
    /// group differs by experimental design: (condition, replicate) for label-free, (file, channel)
    /// for isobaric, and the source file when there is no design at all.
    /// </summary>
    /// <remarks>
    /// Matching a sample group across the records of a dataset must key on this, never on
    /// <see cref="Label"/> (not unique) nor on position within a record's list (a record only has
    /// sample groups for the files it was actually observed in, so the same index means different
    /// files for different records — which silently files one record's counts under another
    /// record's column).
    /// </remarks>
    public string Identity { get; init; } = string.Empty;

    /// <summary>
    /// The file path <see cref="Label"/> was derived from, or null when the label came from the
    /// experimental design instead. Used to widen a label with parent directories when two sample
    /// groups in a dataset would otherwise present the same column name.
    /// </summary>
    public string? LabelSourcePath { get; init; }

    /// <summary>
    /// Identity of the section this result's <see cref="SpectralCount"/> and count-based occupancy
    /// belong to, which is not always the section its intensity belongs to.
    ///
    /// A spectral count answers "how many spectra were acquired", so it is a property of the
    /// acquired FILE. An intensity answers "how much was in this sample", so for an isobaric
    /// experiment it is a property of the CHANNEL, and one file carries many. Reporting a count per
    /// channel restates the file's count once per channel: an 11-plex protein row carries eleven
    /// identical <c>SpectralCount_</c> values and eleven byte-identical <c>CountOccupancy_</c>
    /// strings, three of which sit beside a blank intensity while claiming <c>fraction=1.00(1/1)</c>.
    ///
    /// So the two spaces are declared separately and the schema is told which is which, rather than
    /// inferring one from the other. Defaults to <see cref="Identity"/>, which is correct wherever
    /// the two coincide: a label-free sample group and a design-less file are each their own count
    /// section, and neither their columns nor their values move.
    /// </summary>
    public string CountIdentity
    {
        get => _countIdentity ?? Identity;
        init => _countIdentity = value;
    }

    /// <summary>
    /// Display label for the count section, defaulting to <see cref="Label"/>. Like
    /// <see cref="Label"/> it is not unique and is disambiguated before it reaches a column name,
    /// widened with <see cref="LabelSourcePath"/> -- there is no separate path for the count
    /// section, because the section a result counts under is always a file the result came from.
    /// </summary>
    public string CountLabel
    {
        get => _countLabel ?? Label;
        init => _countLabel = value;
    }

    private readonly string? _countIdentity;
    private readonly string? _countLabel;

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
    /// Formats occupancy for a TSV cell.
    /// Output: semicolon-separated mod entries within each entity, pipe-separated between entities.
    /// </summary>
    /// <param name="orderedKeys">Ordered accessions (protein-level) or base sequences (peptide-level).</param>
    /// <param name="proteinLevel">True for protein-level occupancy; false for peptide-level.</param>
    /// <param name="intensityBased">True to format intensity-based stoichiometry; false for count-based occupancy.</param>
    public string FormatOccupancy(IEnumerable<string> orderedKeys, bool proteinLevel = true, bool intensityBased = false)
    {
        var occupancy = proteinLevel ? ParentOccupancy : DigestionProductOccupancy;

        // A site with no measured intensity has no stoichiometry to report. Formatting it anyway
        // prints "fraction=0.0000(0/0)", which reads as a measured zero rather than as absent data —
        // and does so in rows whose Intensity cell is blank, so the two disagree. Count-based
        // occupancy needs no such filter: its denominator is the observation count, always real.
        Func<SiteSpecificModificationOccupancy, bool> hasEvidence = intensityBased
            ? site => site.TotalIntensity > 0
            : _ => true;

        return FormatOccupancy(occupancy, orderedKeys, hasEvidence, o => o.ToModInfoString(intensityBased));
    }

    /// <summary>
    /// Core formatting helper. Iterates ordered keys, formats each entity's modifications,
    /// and joins with the standard separators (; within entity, | between entities).
    /// </summary>
    private static string FormatOccupancy(
        Dictionary<string, Dictionary<int, List<SiteSpecificModificationOccupancy>>> occupancy,
        IEnumerable<string> orderedKeys,
        Func<SiteSpecificModificationOccupancy, bool> include,
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
                    .Where(include)
                    .Select(formatter));

            if (!string.IsNullOrEmpty(entityString))
                parts.Add(entityString);
        }

        return string.Join("|", parts);
    }

    #endregion
}