using Omics.Fragmentation;

namespace Omics.SpectralMatch
{
    /// <summary>
    /// A spectral match implementation that supports fragment-level sequence coverage calculation.
    /// Implements both <see cref="ISpectralMatch"/> for basic PSM functionality and 
    /// <see cref="IHasSequenceCoverageFromFragments"/> for determining which residues are 
    /// covered by matched fragment ions.
    /// 
    /// Coverage is determined by analyzing N-terminal and C-terminal fragment positions to identify
    /// residues that have fragment evidence on both sides. This works for any biopolymer type
    /// (peptides use b/y ions, nucleic acids use 5'/3' fragments).
    /// </summary>
    public abstract class BaseSpectralMatch : ISpectralMatch, IHasSequenceCoverageFromFragments, IEquatable<BaseSpectralMatch>
    {
        public const double ToleranceForScoreDifferentiation = 1e-9;

        /// <summary>
        /// Creates a new spectral match with the specified properties.
        /// </summary>
        /// <param name="fullFilePath">The file path or identifier for the source spectra file. Null values are converted to empty string.</param>
        /// <param name="oneBasedScanNumber">The one-based scan number for this identification.</param>
        /// <param name="score">The numeric score for this match (higher is better by convention).</param>
        /// <param name="fullSequence">The full modified sequence string with modification annotations. Null values are converted to empty string.</param>
        /// <param name="baseSequence">The unmodified base sequence. Null values are converted to empty string.</param>
        /// <param name="matchedIons"></param>
        protected BaseSpectralMatch(
            string fullFilePath,
            int oneBasedScanNumber,
            double score,
            string fullSequence,
            string baseSequence,
            List<MatchedFragmentIon>? matchedIons = null)
        {
            FullFilePath = fullFilePath ?? string.Empty;
            OneBasedScanNumber = oneBasedScanNumber;
            Score = score;
            FullSequence = fullSequence ?? string.Empty;
            BaseSequence = baseSequence ?? string.Empty;

            MatchedFragmentIons = matchedIons ?? new List<MatchedFragmentIon>();
        }

        #region ISpectralMatch Implementation

        /// <summary>
        /// The file path or identifier for the source spectra file.
        /// Used to associate this match with its originating data file.
        /// </summary>
        /// <inheritdoc cref="ISpectralMatch.FullFilePath"/>
        public string FullFilePath { get; protected set; }

        /// <inheritdoc cref="ISpectralMatch.IsDecoy"/>
        public bool IsDecoy { get; protected set; }

        /// <inheritdoc cref="ISpectralMatch.Accession"/>
        public string Accession { get; protected set; }

        /// <summary>
        /// The full modified sequence string with modification annotations.
        /// For peptides, this includes bracket notation for modifications (e.g., "PEP[Phospho]TIDE").
        /// For nucleic acids, this includes modification annotations in the appropriate format.
        /// </summary>
        /// <inheritdoc cref="ISpectralMatch.FullSequence"/>
        public string FullSequence { get; protected set; }

        /// <summary>
        /// The unmodified base sequence containing only the residue letters.
        /// For peptides, this is the amino acid sequence without any modification annotations.
        /// For nucleic acids, this is the nucleotide sequence.
        /// </summary>
        /// <inheritdoc cref="ISpectralMatch.BaseSequence"/>
        public string BaseSequence { get; protected set; }

        /// <summary>
        /// The one-based scan number for this identification.
        /// This corresponds to the scan index in the source spectra file (first scan = 1).
        /// </summary>
        /// <inheritdoc cref="ISpectralMatch.OneBasedScanNumber"/>
        public int OneBasedScanNumber { get; protected set; }

        /// <summary>
        /// The numeric score for this match. Higher values indicate better matches by convention.
        /// The specific scoring algorithm and scale depend on the search engine used.
        /// </summary>
        /// <inheritdoc cref="ISpectralMatch.Score"/>
        public double Score { get; protected set; }


        /// <inheritdoc cref="ISpectralMatch.Intensities"/>
        public double[]? Intensities { get; protected set; }

        /// <summary>
        /// Returns the biopolymers (peptides, oligonucleotides, etc.) identified for this spectral match.
        /// For ambiguous matches, this may return multiple candidates.
        /// The returned enumerable is a snapshot of the current identifications.
        /// </summary>
        /// <returns>An enumerable of identified biopolymers with their modifications.</returns>
        /// <inheritdoc cref="ISpectralMatch.GetIdentifiedBioPolymersWithSetMods"/>
        public abstract IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods();

        #endregion

        #region IHasSequenceCoverageFromFragments Properties

        public List<MatchedFragmentIon> MatchedFragmentIons { get; set; }

        /// <summary>
        /// One-based positions of residues covered by matched fragment ions.
        /// Null until <see cref="GetSequenceCoverage(IEnumerable{int}, IEnumerable{int})"/> is called.
        /// A residue is considered covered if it has fragment evidence from sequential ions
        /// or from both N-terminal and C-terminal directions.
        /// </summary>
        /// <inheritdoc cref="IHasSequenceCoverageFromFragments.FragmentCoveragePositionInPeptide"/>
        public List<int> FragmentCoveragePositionInPeptide { get; set; }

        public void GetSequenceCoverage()
        {
            FragmentCoveragePositionInPeptide = HasSequenceCoverageFromFragmentsExtensions.GetSequenceCoverage(this);
        }

        #endregion

        #region IComparable Implementation

        /// <summary>
        /// Compares this spectral match to another for ordering purposes.
        /// Orders by Score (ascending), then OneBasedScanNumber (descending).
        /// This ordering ensures higher-scoring matches appear last, with ties broken deterministically.
        /// </summary>
        /// <param name="other">The other spectral match to compare to.</param>
        /// <returns>Negative if this should sort before other; positive if after; zero if equal.</returns>
        public int CompareTo(ISpectralMatch? other)
        {
            if (other is null) return 1;

            if (Math.Abs(this.Score - other.Score) > ToleranceForScoreDifferentiation)
                return this.Score.CompareTo(other.Score);

            // MetaMorpheus would check delta score and precursor mass error here. Since those properties are not part of the interface, we skip them.

            return other.OneBasedScanNumber.CompareTo(OneBasedScanNumber); //reverse the comparision so that the lower scan number comes first.
        }

        /// <summary>
        /// Compares to another <see cref="IHasSequenceCoverageFromFragments"/>.
        /// If the other object also implements <see cref="ISpectralMatch"/>, delegates to 
        /// <see cref="CompareTo(ISpectralMatch)"/>. Otherwise returns 0 (equal) or -1 if other is null.
        /// </summary>
        /// <param name="other">The other object to compare to.</param>
        /// <returns>Comparison result for ordering.</returns>
        public int CompareTo(IHasSequenceCoverageFromFragments? other)
        {
            if (other is ISpectralMatch spectralMatch)
            {
                return CompareTo(spectralMatch);
            }
            return other is null ? -1 : 0;
        }

        #endregion

        #region Equality Implementation

        /// <summary>
        /// Determines equality based on FullFilePath, OneBasedScanNumber, and FullSequence.
        /// Score is intentionally excluded from equality comparison to allow different scoring
        /// methods to identify the same underlying match.
        /// </summary>
        /// <param name="other">The other spectral match to compare.</param>
        /// <returns>True if the matches represent the same identification; otherwise false.</returns>
        public bool Equals(BaseSpectralMatch? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;

            return string.Equals(FullFilePath, other.FullFilePath, StringComparison.Ordinal)
                && OneBasedScanNumber == other.OneBasedScanNumber
                && string.Equals(FullSequence, other.FullSequence, StringComparison.Ordinal);
        }

        /// <summary>
        /// Determines equality with another object.
        /// Only returns true for <see cref="BaseSpectralMatch"/> instances
        /// that match on FullFilePath, OneBasedScanNumber, and FullSequence.
        /// </summary>
        /// <param name="obj">The object to compare.</param>
        /// <returns>True if equal; otherwise false.</returns>
        public override bool Equals(object? obj)
        {
            if (obj is BaseSpectralMatch sm) return Equals(sm);
            return false;
        }

        /// <summary>
        /// Returns a hash code based on FullFilePath, OneBasedScanNumber, and FullSequence.
        /// Consistent with <see cref="Equals(BaseSpectralMatch)"/> implementation.
        /// </summary>
        /// <returns>A hash code for this instance.</returns>
        public override int GetHashCode()
        {
            return HashCode.Combine(FullFilePath, OneBasedScanNumber, FullSequence);
        }

        /// <summary>
        /// Equality operator. Two matches are equal if they have the same FullFilePath, 
        /// OneBasedScanNumber, and FullSequence.
        /// </summary>
        public static bool operator ==(BaseSpectralMatch? left, BaseSpectralMatch? right)
        {
            if (left is null) return right is null;
            return left.Equals(right);
        }

        /// <summary>
        /// Inequality operator.
        /// </summary>
        public static bool operator !=(BaseSpectralMatch? left, BaseSpectralMatch? right)
        {
            return !(left == right);
        }

        #endregion

        #region Object Overrides

        /// <summary>
        /// Returns a human-readable string representation of this spectral match.
        /// Format: "Scan {number}: {sequence} (Score: {score:F2})"
        /// </summary>
        /// <returns>A string representation suitable for debugging and logging.</returns>
        public override string ToString()
        {
            return $"Scan {OneBasedScanNumber}: {FullSequence} (Score: {Score:F2})";
        }

        #endregion
    }
}