using MassSpectrometry;
using Omics.Modifications;

namespace Omics.BioPolymerGroup
{
    /// <summary>
    /// Interface for biopolymer groups that supports both label-free and isobaric quantification.
    /// A biopolymer group represents a collection of related biopolymers (e.g., proteins, oligonucleotides)
    /// that share peptides/fragments and are grouped together for quantification and statistical analysis.
    /// </summary>
    public interface IBioPolymerGroup : IEquatable<IBioPolymerGroup>
    {
        /// <summary>
        /// True if this group contains only decoy biopolymers, used for FDR estimation.
        /// </summary>
        bool IsDecoy { get; }

        /// <summary>
        /// True if this group contains biopolymers marked as contaminants.
        /// </summary>
        bool IsContaminant { get; }

        /// <summary>
        /// List of files or samples that contribute quantification data for this group.
        /// Supports both SpectraFileInfo (label-free) and IsobaricQuantSampleInfo (TMT/iTRAQ).
        /// </summary>
        List<ISampleInfo> SamplesForQuantification { get; set; }

        /// <summary>
        /// Set of all biopolymers (e.g., proteins, RNA sequences) that belong to this group.
        /// </summary>
        HashSet<IBioPolymer> BioPolymers { get; set; }

        /// <summary>
        /// Display name for the biopolymer group, typically derived from accession(s) or gene name(s).
        /// Used as the primary identity key for equality comparisons.
        /// </summary>
        string BioPolymerGroupName { get; }

        /// <summary>
        /// Aggregated score for the group, computed from member PSMs or peptides.
        /// Higher scores typically indicate higher confidence.
        /// </summary>
        double BioPolymerGroupScore { get; set; }

        /// <summary>
        /// All biopolymer sequences with set modifications identified in this group,
        /// including those shared with other groups.
        /// </summary>
        HashSet<IBioPolymerWithSetMods> AllBioPolymersWithSetMods { get; set; }

        /// <summary>
        /// Biopolymer sequences with set modifications that are unique to this group
        /// (not shared with any other biopolymer group).
        /// </summary>
        HashSet<IBioPolymerWithSetMods> UniqueBioPolymersWithSetMods { get; set; }

        /// <summary>
        /// All peptide-spectrum matches (PSMs) for this group that pass the 1% FDR threshold.
        /// </summary>
        HashSet<ISpectralMatch> AllPsmsBelowOnePercentFDR { get; set; }

        /// <summary>
        /// The q-value (FDR-adjusted p-value) for this biopolymer group.
        /// Lower values indicate higher confidence in the identification.
        /// </summary>
        double QValue { get; set; }

        /// <summary>
        /// The best (lowest) q-value among all biopolymers with set modifications in this group.
        /// </summary>
        double BestBioPolymerWithSetModsQValue { get; set; }

        /// <summary>
        /// The best (highest) score among all biopolymers with set modifications in this group.
        /// </summary>
        double BestBioPolymerWithSetModsScore { get; set; }

        /// <summary>
        /// Summary information about modifications present on members of this group.
        /// Each string typically describes a modification type and its frequency or location.
        /// </summary>
        List<string> ModsInfo { get; }

        /// <summary>
        /// Dictionary mapping sample identifiers to measured intensity values for this group.
        /// Supports both SpectraFileInfo and IsobaricQuantSampleInfo as keys.
        /// </summary>
        Dictionary<ISampleInfo, double> IntensitiesBySample { get; set; }

        /// <summary>
        /// All biopolymers in this group ordered alphabetically by accession.
        /// Provides a stable, deterministic ordering for output and comparison.
        /// </summary>
        List<IBioPolymer> ListOfBioPolymersOrderedByAccession { get; }

        /// <summary>
        /// Returns a tab-separated header line for output files, matching the format of <see cref="ToString"/>.
        /// </summary>
        string GetTabSeparatedHeader();

        /// <summary>
        /// Computes or updates the <see cref="BioPolymerGroupScore"/> based on member PSMs and peptides.
        /// </summary>
        void Score();

        /// <summary>
        /// Merges another biopolymer group into this one, combining their members, PSMs, and intensities.
        /// Used when groups are determined to represent the same biological entity.
        /// </summary>
        /// <param name="otherBioPolymerGroup">The group to merge into this one.</param>
        void MergeWith(IBioPolymerGroup otherBioPolymerGroup);

        /// <summary>
        /// Creates a new biopolymer group containing only data from a specific file.
        /// Used for per-file analysis and output.
        /// </summary>
        /// <param name="fullFilePath">The full path to the file to subset by.</param>
        /// <param name="silacLabels">Optional SILAC labels to apply during subset creation.</param>
        /// <returns>A new group containing only data from the specified file.</returns>
        IBioPolymerGroup ConstructSubsetBioPolymerGroup(string fullFilePath, List<SilacLabel>? silacLabels = null);
    }
}