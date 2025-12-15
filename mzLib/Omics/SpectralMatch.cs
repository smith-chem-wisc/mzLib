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
        /// Positions in the biopolymer sequence (one-based) that are covered by fragment ions.
        /// Populated by <see cref="GetSequenceCoverage"/> when fragment coverage data is available.
        /// May be null if coverage has not been calculated or is not available.
        /// </summary>
        public HashSet<int>? FragmentCoveragePositionInPeptide { get; protected set; }

        /// <summary>
        /// N-terminal (or 5' for nucleic acids) fragment positions (one-based).
        /// Derived classes should populate this before calling <see cref="GetSequenceCoverage"/>.
        /// </summary>
        protected List<int>? NTerminalFragmentPositions { get; set; }

        /// <summary>
        /// C-terminal (or 3' for nucleic acids) fragment positions (one-based).
        /// Derived classes should populate this before calling <see cref="GetSequenceCoverage"/>.
        /// </summary>
        protected List<int>? CTerminalFragmentPositions { get; set; }

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
        /// Calculates sequence coverage from fragment ions for this spectral match.
        /// Populates <see cref="FragmentCoveragePositionInPeptide"/> with one-based positions
        /// of residues that are covered by matched fragment ions.
        /// Works for any biopolymer type (proteins, nucleic acids, etc.).
        /// 
        /// Derived classes should populate <see cref="NTerminalFragmentPositions"/> and 
        /// <see cref="CTerminalFragmentPositions"/> before calling this method, or override
        /// this method entirely to provide custom coverage calculation.
        /// </summary>
        public virtual void GetSequenceCoverage()
        {
            if (string.IsNullOrEmpty(BaseSequence))
            {
                return;
            }

            var nTermPositions = NTerminalFragmentPositions ?? new List<int>();
            var cTermPositions = CTerminalFragmentPositions ?? new List<int>();

            if (!nTermPositions.Any() && !cTermPositions.Any())
            {
                return;
            }

            var fragmentCoveredResidues = new HashSet<int>();

            // Process N-terminal fragments (or 5' for nucleic acids)
            if (nTermPositions.Any())
            {
                var sortedNTerm = nTermPositions.OrderBy(x => x).ToList();

                // If the final N-terminal fragment is present, last residue is covered
                if (sortedNTerm.Contains(BaseSequence.Length - 1))
                {
                    fragmentCoveredResidues.Add(BaseSequence.Length);
                }

                // If the first N-terminal fragment is present, first residue is covered
                if (sortedNTerm.Contains(1))
                {
                    fragmentCoveredResidues.Add(1);
                }

                // Check all positions except for the last one in the list
                for (int i = 0; i < sortedNTerm.Count - 1; i++)
                {
                    // Sequential positions mean the second one is covered
                    if (sortedNTerm[i + 1] - sortedNTerm[i] == 1)
                    {
                        fragmentCoveredResidues.Add(sortedNTerm[i + 1]);
                    }

                    // Check if position is covered from both directions (inclusive)
                    if (cTermPositions.Contains(sortedNTerm[i + 1]))
                    {
                        fragmentCoveredResidues.Add(sortedNTerm[i + 1]);
                    }

                    // Check if position is covered from both directions (exclusive)
                    if (cTermPositions.Contains(sortedNTerm[i + 1] + 2))
                    {
                        fragmentCoveredResidues.Add(sortedNTerm[i + 1] + 1);
                    }
                }
            }

            // Process C-terminal fragments (or 3' for nucleic acids)
            if (cTermPositions.Any())
            {
                var sortedCTerm = cTermPositions.OrderBy(x => x).ToList();

                // If the second position is present, first residue is covered
                if (sortedCTerm.Contains(2))
                {
                    fragmentCoveredResidues.Add(1);
                }

                // If the last position is present, final residue is covered
                if (sortedCTerm.Contains(BaseSequence.Length))
                {
                    fragmentCoveredResidues.Add(BaseSequence.Length);
                }

                // Check all positions except for the last one in the list
                for (int i = 0; i < sortedCTerm.Count - 1; i++)
                {
                    // Sequential positions mean the first one is covered
                    if (sortedCTerm[i + 1] - sortedCTerm[i] == 1)
                    {
                        fragmentCoveredResidues.Add(sortedCTerm[i]);
                    }
                }
            }

            FragmentCoveragePositionInPeptide = fragmentCoveredResidues;
        }

        /// <summary>
        /// Sets the fragment positions for coverage calculation.
        /// Call this method before <see cref="GetSequenceCoverage"/> to provide fragment data.
        /// </summary>
        /// <param name="nTerminalPositions">One-based positions of N-terminal (or 5') fragments.</param>
        /// <param name="cTerminalPositions">One-based positions of C-terminal (or 3') fragments.</param>
        public void SetFragmentPositions(IEnumerable<int>? nTerminalPositions, IEnumerable<int>? cTerminalPositions)
        {
            NTerminalFragmentPositions = nTerminalPositions?.ToList();
            CTerminalFragmentPositions = cTerminalPositions?.ToList();
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
        public int CompareTo(ISpectralMatch? other)
        {
            if (other is null) return -1;

            int scoreComparison = other.Score.CompareTo(Score);
            if (scoreComparison != 0) return scoreComparison;

            int fileComparison = string.Compare(FullFilePath, other.FullFilePath, StringComparison.Ordinal);
            if (fileComparison != 0) return fileComparison;

            return OneBasedScanNumber.CompareTo(other.OneBasedScanNumber);
        }

        /// <summary>
        /// Determines whether this spectral match equals another.
        /// Two matches are equal if they have the same FullFilePath, OneBasedScanNumber, and FullSequence.
        /// </summary>
        public bool Equals(SpectralMatch? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;

            return string.Equals(FullFilePath, other.FullFilePath, StringComparison.Ordinal)
                && OneBasedScanNumber == other.OneBasedScanNumber
                && string.Equals(FullSequence, other.FullSequence, StringComparison.Ordinal);
        }

        public override bool Equals(object? obj)
        {
            if (obj is SpectralMatch sm) return Equals(sm);
            return false;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(FullFilePath, OneBasedScanNumber, FullSequence);
        }

        public override string ToString()
        {
            return $"Scan {OneBasedScanNumber}: {FullSequence} (Score: {Score:F2})";
        }

        public static bool operator ==(SpectralMatch? left, SpectralMatch? right)
        {
            if (left is null) return right is null;
            return left.Equals(right);
        }

        public static bool operator !=(SpectralMatch? left, SpectralMatch? right)
        {
            return !(left == right);
        }
    }
}