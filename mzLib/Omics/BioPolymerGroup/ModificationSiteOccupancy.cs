namespace Omics.BioPolymerGroup;

/// <summary>
/// Represents the occupancy/stoichiometry of a single modification at a specific
/// position on a biopolymer. Supports both count-based and intensity-based metrics.
/// </summary>
public class ModificationSiteOccupancy
{
    /// <summary>One-based position in the parent biopolymer sequence.</summary>
    public int OneBasedPositionInBioPolymer { get; }

    /// <summary>The modification identity (e.g., "Oxidation on M").</summary>
    public string ModificationIdWithMotif { get; }

    /// <summary>Number of peptides carrying this mod at this position.</summary>
    public int ModifiedCount { get; set; }

    /// <summary>Total peptides covering this position (modified + unmodified).</summary>
    public int TotalCount { get; set; }

    /// <summary>Count-based occupancy fraction (ModifiedCount / TotalCount).</summary>
    public double CountBasedOccupancy => TotalCount > 0 ? (double)ModifiedCount / TotalCount : 0;

    /// <summary>Sum of intensities for peptides carrying this mod at this position.</summary>
    public double ModifiedIntensity { get; set; }

    /// <summary>Sum of intensities for all peptides covering this position.</summary>
    public double TotalIntensity { get; set; }

    /// <summary>Intensity-based stoichiometry fraction (ModifiedIntensity / TotalIntensity).</summary>
    public double IntensityBasedStoichiometry => TotalIntensity > 0 ? ModifiedIntensity / TotalIntensity : 0;

    public ModificationSiteOccupancy(int oneBasedPosition, string modIdWithMotif)
    {
        OneBasedPositionInBioPolymer = oneBasedPosition;
        ModificationIdWithMotif = modIdWithMotif;
    }

    /// <summary>
    /// Formatted string matching the existing ModsInfo format for backward compatibility.
    /// Format: #aa{position}[{modName},info:occupancy={fraction}({count}/{total})]
    /// </summary>
    public string ToModInfoString()
    {
        string occupancy = CountBasedOccupancy.ToString("F2");
        string fractional = $"{ModifiedCount}/{TotalCount}";
        return $"#aa{OneBasedPositionInBioPolymer}[{ModificationIdWithMotif},info:occupancy={occupancy}({fractional})]";
    }
}
