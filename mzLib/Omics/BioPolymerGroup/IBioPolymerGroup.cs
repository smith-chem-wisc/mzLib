using MassSpectrometry;
using Omics.Modifications;
using Omics.SpectralMatch;

namespace Omics.BioPolymerGroup
{
    /// <summary>
    /// Represents a group of related biopolymers (e.g., proteins, RNA sequences) that cannot be 
    /// distinguished based on the identified peptide or oligonucleotide sequences. Groups are formed 
    /// during protein/gene inference when multiple biopolymers share all their identified sequences.
    /// 
    /// Supports both label-free and isobaric (TMT/iTRAQ) quantification methods.
    /// Implementations provide scoring, FDR estimation, sequence coverage, and output formatting.
    /// </summary>
    public interface IBioPolymerGroup : IEquatable<IBioPolymerGroup>
    {
        /// <summary>
        /// True if any biopolymer in this group is marked as a decoy, used for FDR estimation.
        /// </summary>
        bool IsDecoy { get; }

        /// <summary>
        /// True if any biopolymer in this group is marked as a contaminant.
        /// </summary>
        bool IsContaminant { get; }

        /// <summary>
        /// Samples that contribute quantification data for this group.
        /// Supports <see cref="SpectraFileInfo"/> (label-free) and <see cref="IsobaricQuantSampleInfo"/> (TMT/iTRAQ).
        /// </summary>
        List<ISampleInfo> SamplesForQuantification { get; set; }

        /// <summary>
        /// All biopolymers (e.g., proteins, RNA sequences) that belong to this group.
        /// These biopolymers are indistinguishable based on the identified sequences.
        /// </summary>
        HashSet<IBioPolymer> BioPolymers { get; set; }

        /// <summary>
        /// Display name for the group, typically a pipe-delimited concatenation of member accessions.
        /// Used as the primary identity key for equality comparisons via <see cref="IEquatable{T}"/>.
        /// </summary>
        string BioPolymerGroupName { get; }

        /// <summary>
        /// Aggregated confidence score for the group. Higher values indicate higher confidence.
        /// 
        /// Computed by <see cref="BioPolymerGroupExtensions.Score"/> as the sum of the best (highest) 
        /// score for each unique base sequence among the PSMs in <see cref="AllPsmsBelowOnePercentFDR"/>.
        /// This ensures each unique peptide/oligonucleotide contributes only its best-scoring identification.
        /// </summary>
        /// <seealso cref="BioPolymerGroupExtensions.Score"/>
        double BioPolymerGroupScore { get; set; }

        /// <summary>
        /// All identified sequences with modifications in this group, including sequences 
        /// that are shared with other biopolymer groups.
        /// </summary>
        HashSet<IBioPolymerWithSetMods> AllBioPolymersWithSetMods { get; set; }

        /// <summary>
        /// Identified sequences with modifications that are unique to this group 
        /// (not found in any other biopolymer group). Used for protein inference.
        /// </summary>
        HashSet<IBioPolymerWithSetMods> UniqueBioPolymersWithSetMods { get; set; }

        /// <summary>
        /// Peptide-spectrum matches (PSMs) for this group that pass the 1% FDR threshold.
        /// Used for scoring, coverage calculation, and quantification.
        /// </summary>
        HashSet<ISpectralMatch> AllPsmsBelowOnePercentFDR { get; set; }

        /// <summary>
        /// The q-value for this biopolymer group, representing the minimum FDR at which 
        /// this group would be accepted. Lower values indicate higher confidence (0.01 = 1% FDR).
        /// </summary>
        double QValue { get; set; }

        /// <summary>
        /// The best (lowest) q-value among all identified sequences in this group.
        /// </summary>
        double BestBioPolymerWithSetModsQValue { get; set; }

        /// <summary>
        /// The best (highest) score among all identified sequences in this group.
        /// </summary>
        double BestBioPolymerWithSetModsScore { get; set; }

        /// <summary>
        /// Measured intensity values for this group, keyed by sample.
        /// Supports both <see cref="SpectraFileInfo"/> and <see cref="IsobaricQuantSampleInfo"/> as keys.
        /// </summary>
        Dictionary<ISampleInfo, double> IntensitiesBySample { get; set; }

        /// <summary>
        /// All biopolymers in this group ordered alphabetically by accession.
        /// Provides stable, deterministic ordering for output and display.
        /// </summary>
        List<IBioPolymer> ListOfBioPolymersOrderedByAccession { get; }

        /// <summary>
        /// Returns a tab-separated header line for output files.
        /// The format matches the output of <see cref="object.ToString"/>.
        /// </summary>
        /// <returns>Tab-separated header string suitable for TSV file output.</returns>
        string GetTabSeparatedHeader();

        /// <summary>
        /// Merges another biopolymer group into this one, combining members, PSMs, sequences, 
        /// and intensities. The merged group's score is reset to 0.
        /// Used when groups are determined to represent the same biological entity.
        /// </summary>
        /// <param name="otherBioPolymerGroup">The group to merge into this one.</param>
        void MergeWith(IBioPolymerGroup otherBioPolymerGroup);

        /// <summary>
        /// Creates a new biopolymer group containing only data from a specific spectra file.
        /// The new group shares the same biopolymers but has filtered PSMs, sequences, and intensities.
        /// </summary>
        /// <param name="fullFilePath">The full path to the spectra file to filter by.</param>
        /// <param name="silacLabels">Optional SILAC labels to apply during subset creation.</param>
        /// <returns>A new <see cref="IBioPolymerGroup"/> containing only data from the specified file.</returns>
        IBioPolymerGroup ConstructSubsetBioPolymerGroup(string fullFilePath, List<SilacLabel>? silacLabels = null);
    }

    /// <summary>
    /// Extension methods for <see cref="IBioPolymerGroup"/> that provide common operations
    /// applicable to all implementations.
    /// </summary>
    public static class BioPolymerGroupExtensions
    {
        /// <summary>
        /// Calculates and updates the <see cref="IBioPolymerGroup.BioPolymerGroupScore"/> based on PSM scores.
        /// 
        /// The score is computed as the sum of the best (highest) score for each unique base sequence
        /// among the PSMs in <see cref="IBioPolymerGroup.AllPsmsBelowOnePercentFDR"/>. This approach
        /// ensures that each unique peptide/oligonucleotide sequence contributes only its best-scoring
        /// identification to the group score.
        /// </summary>
        /// <param name="group">The biopolymer group to score.</param>
        /// <remarks>
        /// This method should be called after <see cref="IBioPolymerGroup.AllPsmsBelowOnePercentFDR"/> 
        /// has been populated. If the collection is empty, the score will be set to 0.
        /// </remarks>
        /// <example>
        /// <code>
        /// IBioPolymerGroup group = GetBioPolymerGroup();
        /// group.Score(); // Calculates and sets BioPolymerGroupScore
        /// </code>
        /// </example>
        public static void Score(this IBioPolymerGroup group)
        {
            group.BioPolymerGroupScore = group.AllPsmsBelowOnePercentFDR
                .GroupBy(p => p.BaseSequence)
                .Select(p => p.Select(x => x.Score).Max())
                .Sum();
        }
    }
}