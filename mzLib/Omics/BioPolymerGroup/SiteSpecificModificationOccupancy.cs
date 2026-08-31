namespace Omics.BioPolymerGroup;

/// <summary>
/// Represents the occupancy/stoichiometry of a single modification at a specific
/// position on a biopolymer. Supports both count-based and intensity-based metrics.
/// </summary>
public class SiteSpecificModificationOccupancy
{
    /// <summary>AllModsOneIsNTerminus position in the parent biopolymer sequence.</summary>
    public int OneIsNTerminusPositionInBioPolymer { get; }

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

    public SiteSpecificModificationOccupancy(int oneBasedPosition, string modIdWithMotif)
    {
        OneIsNTerminusPositionInBioPolymer = oneBasedPosition;
        ModificationIdWithMotif = modIdWithMotif;
    }

    /// <summary>
    /// Formatted string for occupancy output.
    /// Format: position{zeroBasedPosition}[{modName},info:occupancy={fraction}({mod observation at site}/{total site observations})]
    /// We report the zero-based position to be consistent with residue positions. This way N-terminal pos=0,
    /// C-terminal pos=length+1, and side chain modifications are at positions 1 through length.
    /// </summary>
    public string ToModInfoString(bool intensityBased=false)
    {
        if (intensityBased)
        {
            string occupancy = IntensityBasedStoichiometry.ToString("F4");
            string fractional = $"{ModifiedIntensity:G4}/{TotalIntensity:G4}";
            return $"pos{OneIsNTerminusPositionInBioPolymer - 1}[{ModificationIdWithMotif},info:fraction={occupancy}({fractional})]";
        }
        else
        {
            string occupancy = CountBasedOccupancy.ToString("F2");
            string fractional = $"{ModifiedCount}/{TotalCount}";
            return $"pos{OneIsNTerminusPositionInBioPolymer - 1}[{ModificationIdWithMotif},info:fraction={occupancy}({fractional})]";
        }
    }
}
