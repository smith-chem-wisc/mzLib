namespace Omics
{
    /// <summary>
    /// Represents a single spectrum-to-sequence identification (PSM).
    /// This concrete implementation of <see cref="ISpectralMatch"/> provides basic storage
    /// and comparison for spectral match data produced by search engines or readers.
    /// </summary>
    public class SpectralMatch : ISpectralMatch, IEquatable<SpectralMatch>
    {
        /// <summary>
        /// Creates a new spectral match with the specified properties.
        /// </summary>
        /// <param name="fullFilePath">The file path or identifier for the source spectra file.</param>
        /// <param name="oneBasedScanNumber">The one-based scan number for this identification.</param>
        /// <param name="score">The numeric score for this match (higher is better).</param>
        /// <param name="fullSequence">The full modified sequence string.</param>
        /// <param name="baseSequence">The unmodified base sequence.</param>
        /// <param name="identifiedBioPolymers">The biopolymers identified for this match.</param>
        public SpectralMatch(
            string fullFilePath,
            int oneBasedScanNumber,
            double score,
            string fullSequence,
            string baseSequence,
            IEnumerable<IBioPolymerWithSetMods>? identifiedBioPolymers = null)
        {
            FullFilePath = fullFilePath ?? string.Empty;
            OneBasedScanNumber = oneBasedScanNumber;
            Score = score;
            FullSequence = fullSequence ?? string.Empty;
            BaseSequence = baseSequence ?? string.Empty;
            _identifiedBioPolymers = identifiedBioPolymers?.ToList() ?? new List<IBioPolymerWithSetMods>();
        }

        /// <summary>
        /// The file path or identifier for the spectra file that produced this match.
        /// </summary>
        public string FullFilePath { get; }

        /// <summary>
        /// The full modified sequence string including modification annotations.
        /// </summary>
        public string FullSequence { get; }

        /// <summary>
        /// The unmodified base sequence.
        /// </summary>
        public string BaseSequence { get; }

        /// <summary>
        /// The one-based scan number for this identification.
        /// </summary>
        public int OneBasedScanNumber { get; }

        /// <summary>
        /// The numeric score for this match. Higher values indicate better matches.
        /// </summary>
        public double Score { get; }

        /// <summary>
        /// Positions in the peptide (one-based) that are covered by fragment ions.
        /// Populated by <see cref="GetAminoAcidCoverage"/> when fragment coverage data is available.
        /// May be null if coverage has not been calculated or is not available.
        /// </summary>
        public HashSet<int>? FragmentCoveragePositionInPeptide { get; protected set; }

        private readonly List<IBioPolymerWithSetMods> _identifiedBioPolymers;

        /// <summary>
        /// Returns the biopolymers identified for this spectral match.
        /// </summary>
        /// <returns>An enumerable of identified biopolymers; never null but may be empty.</returns>
        public IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods()
        {
            return _identifiedBioPolymers;
        }

        /// <summary>
        /// Calculates amino acid coverage from fragment ions for this spectral match.
        /// Derived classes can override this to provide actual fragment coverage calculation.
        /// The base implementation does nothing; override in derived classes that have access to fragment data.
        /// </summary>
        public virtual void GetAminoAcidCoverage()
        {
            // Base implementation does nothing.
            // Derived classes with access to fragment ion data should override this
            // to populate FragmentCoveragePositionInPeptide.
        }

        /// <summary>
        /// Adds a biopolymer to the list of identified biopolymers for this match.
        /// </summary>
        /// <param name="bioPolymer">The biopolymer to add.</param>
        public void AddIdentifiedBioPolymer(IBioPolymerWithSetMods bioPolymer)
        {
            if (bioPolymer != null)
            {
                _identifiedBioPolymers.Add(bioPolymer);
            }
        }

        /// <summary>
        /// Adds multiple biopolymers to the list of identified biopolymers for this match.
        /// </summary>
        /// <param name="bioPolymers">The biopolymers to add.</param>
        public void AddIdentifiedBioPolymers(IEnumerable<IBioPolymerWithSetMods> bioPolymers)
        {
            if (bioPolymers != null)
            {
                _identifiedBioPolymers.AddRange(bioPolymers.Where(b => b != null));
            }
        }

        /// <summary>
        /// Compares this spectral match to another for ordering purposes.
        /// Comparison is by Score (descending), then by FullFilePath, then by OneBasedScanNumber.
        /// </summary>
        /// <param name="other">The other spectral match to compare to.</param>
        /// <returns>
        /// A negative value if this instance precedes <paramref name="other"/>;
        /// zero if they are equal; a positive value if this instance follows <paramref name="other"/>.
        /// </returns>
        public int CompareTo(ISpectralMatch? other)
        {
            if (other is null) return -1;

            // Higher score comes first (descending order)
            int scoreComparison = other.Score.CompareTo(Score);
            if (scoreComparison != 0) return scoreComparison;

            // Then by file path (ascending)
            int fileComparison = string.Compare(FullFilePath, other.FullFilePath, StringComparison.Ordinal);
            if (fileComparison != 0) return fileComparison;

            // Then by scan number (ascending)
            return OneBasedScanNumber.CompareTo(other.OneBasedScanNumber);
        }

        /// <summary>
        /// Determines whether this spectral match equals another.
        /// Two matches are equal if they have the same FullFilePath, OneBasedScanNumber, and FullSequence.
        /// </summary>
        /// <param name="other">The other spectral match to compare.</param>
        /// <returns>True if the matches are equal; otherwise false.</returns>
        public bool Equals(SpectralMatch? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;

            return string.Equals(FullFilePath, other.FullFilePath, StringComparison.Ordinal)
                && OneBasedScanNumber == other.OneBasedScanNumber
                && string.Equals(FullSequence, other.FullSequence, StringComparison.Ordinal);
        }

        /// <summary>
        /// Determines whether this spectral match equals another object.
        /// </summary>
        /// <param name="obj">The object to compare.</param>
        /// <returns>True if the objects are equal; otherwise false.</returns>
        public override bool Equals(object? obj)
        {
            if (obj is SpectralMatch sm) return Equals(sm);
            return false;
        }

        /// <summary>
        /// Returns a hash code based on FullFilePath, OneBasedScanNumber, and FullSequence.
        /// </summary>
        /// <returns>A hash code for this instance.</returns>
        public override int GetHashCode()
        {
            return HashCode.Combine(FullFilePath, OneBasedScanNumber, FullSequence);
        }

        /// <summary>
        /// Returns a string representation of this spectral match.
        /// </summary>
        /// <returns>A string in the format "Scan {ScanNumber}: {FullSequence} (Score: {Score})".</returns>
        public override string ToString()
        {
            return $"Scan {OneBasedScanNumber}: {FullSequence} (Score: {Score:F2})";
        }

        /// <summary>
        /// Determines whether two spectral matches are equal.
        /// </summary>
        public static bool operator ==(SpectralMatch? left, SpectralMatch? right)
        {
            if (left is null) return right is null;
            return left.Equals(right);
        }

        /// <summary>
        /// Determines whether two spectral matches are not equal.
        /// </summary>
        public static bool operator !=(SpectralMatch? left, SpectralMatch? right)
        {
            return !(left == right);
        }
    }
}