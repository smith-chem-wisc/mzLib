using MassSpectrometry;
using Omics.Modifications;

namespace Omics.BioPolymerGroup
{
    /// <summary>
    /// Generic interface for biopolymer groups that allows different quantification strategies.
    /// A biopolymer group represents a collection of related biopolymers (e.g., proteins, oligonucleotides)
    /// that share peptides/fragments and are grouped together for quantification and statistical analysis.
    /// </summary>
    /// <typeparam name="TFileInfo">The file identifier type used to key intensity values
    /// (e.g., <see cref="SpectraFileInfo"/> for label-free, <see cref="IsobaricQuantFileInfo"/> for isobaric).</typeparam>
    /// <typeparam name="TIntensity">The intensity value type
    /// (e.g., <see cref="double"/> for label-free single values, <see cref="double"/>[] for isobaric channel arrays).</typeparam>
    internal interface IBioPolymerGroup<TFileInfo, TIntensity> : IEquatable<IBioPolymerGroup<TFileInfo, TIntensity>>
        where TFileInfo : notnull
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
        /// List of files that contribute quantification data for this group.
        /// The type depends on the quantification strategy (label-free vs. isobaric).
        /// </summary>
        List<TFileInfo> FilesForQuantification { get; set; }

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
        /// All unique biopolymer sequences with set modifications identified in this group,
        /// including those shared with other groups.
        /// </summary>
        HashSet<IBioPolymerWithSetMods> AllBioPolymerWithSetMods { get; set; }

        /// <summary>
        /// Biopolymer sequences with set modifications that are unique to this group
        /// (not shared with any other biopolymer group).
        /// </summary>
        HashSet<IBioPolymerWithSetMods> UniqueBioPolymerWithSetMods { get; set; }

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
        double BestBiopolymerWithSetsModQValue { get; set; }

        /// <summary>
        /// The best (highest) score among all biopolymers with set modifications in this group.
        /// </summary>
        double BestBioPolymerWithSetsModScore { get; set; }

        /// <summary>
        /// Summary information about modifications present on members of this group.
        /// Each string typically describes a modification type and its frequency or location.
        /// </summary>
        List<string> ModsInfo { get; }

        /// <summary>
        /// Dictionary mapping file identifiers to measured intensity values for this group.
        /// For label-free quantification, values are single doubles.
        /// For isobaric quantification (TMT/iTRAQ), values are arrays of channel intensities.
        /// </summary>
        Dictionary<TFileInfo, TIntensity> IntensitiesByFile { get; set; }

        /// <summary>
        /// All biopolymers in this group ordered alphabetically by accession.
        /// Provides a stable, deterministic ordering for output and comparison.
        /// </summary>
        List<IBioPolymer> ListOfBioPolymersOrderedByAccession { get; }

        /// <summary>
        /// Default equality implementation based on <see cref="BioPolymerGroupName"/>.
        /// Two groups are equal if they have the same name.
        /// </summary>
        bool IEquatable<IBioPolymerGroup<TFileInfo, TIntensity>>.Equals(IBioPolymerGroup<TFileInfo, TIntensity>? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;
            if (other.GetType() != GetType()) return false;
            return BioPolymerGroupName == other.BioPolymerGroupName;
        }

        /// <summary>
        /// Returns a tab-separated header line for output files, matching the format of <see cref="ToString"/>.
        /// </summary>
        string GetTabSeparatedHeader();

        /// <summary>
        /// Returns a tab-separated string representation of this group for output files.
        /// </summary>
        string ToString();

        /// <summary>
        /// Computes or updates the <see cref="BioPolymerGroupScore"/> based on member PSMs and peptides.
        /// </summary>
        void Score();

        /// <summary>
        /// Merges another biopolymer group into this one, combining their members, PSMs, and intensities.
        /// Used when groups are determined to represent the same biological entity.
        /// </summary>
        /// <param name="otherBioPolymerGroup">The group to merge into this one.</param>
        void MergeWith(IBioPolymerGroup<TFileInfo, TIntensity> otherBioPolymerGroup);

        /// <summary>
        /// Creates a new biopolymer group containing only data from a specific file.
        /// Used for per-file analysis and output.
        /// </summary>
        /// <param name="fullFilePath">The full path to the file to subset by.</param>
        /// <param name="silacLabels">Optional SILAC labels to apply during subsetting.</param>
        /// <returns>A new group containing only data from the specified file.</returns>
        IBioPolymerGroup<TFileInfo, TIntensity> ConstructSubsetBioPolymerGroup(string fullFilePath, List<SilacLabel> silacLabels = null);
    }

    /// <summary>
    /// Convenience interface for label-free quantification biopolymer groups.
    /// Uses <see cref="SpectraFileInfo"/> as the file key and <see cref="double"/> for single intensity values.
    /// </summary>
    internal interface ILabelFreeBioPolymerGroup : IBioPolymerGroup<SpectraFileInfo, double> { }

    /// <summary>
    /// Convenience interface for isobaric (TMT/iTRAQ) quantification biopolymer groups.
    /// Uses <see cref="IsobaricQuantFileInfo"/> as the file key and <see cref="double"/>[] for per-channel intensities.
    /// </summary>
    internal interface IIsobaricBioPolymerGroup : IBioPolymerGroup<IsobaricQuantFileInfo, double[]> { }
}