namespace Omics
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
    public class SpectralMatchWithSequenceCoverage : ISpectralMatch, IHasSequenceCoverageFromFragments, IEquatable<SpectralMatchWithSequenceCoverage>
    {
        /// <summary>
        /// Creates a new spectral match with the specified properties.
        /// </summary>
        /// <param name="fullFilePath">The file path or identifier for the source spectra file. Null values are converted to empty string.</param>
        /// <param name="oneBasedScanNumber">The one-based scan number for this identification.</param>
        /// <param name="score">The numeric score for this match (higher is better by convention).</param>
        /// <param name="fullSequence">The full modified sequence string with modification annotations. Null values are converted to empty string.</param>
        /// <param name="baseSequence">The unmodified base sequence. Null values are converted to empty string.</param>
        /// <param name="identifiedBioPolymers">Optional biopolymers identified for this match. A defensive copy is made.</param>
        public SpectralMatchWithSequenceCoverage(
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

        /// <inheritdoc cref="ISpectralMatch.FullFilePath"/>
        public string FullFilePath { get; }

        /// <inheritdoc cref="ISpectralMatch.FullSequence"/>
        public string FullSequence { get; }

        /// <inheritdoc cref="ISpectralMatch.BaseSequence"/>
        public string BaseSequence { get; }

        /// <inheritdoc cref="ISpectralMatch.OneBasedScanNumber"/>
        public int OneBasedScanNumber { get; }

        /// <inheritdoc cref="ISpectralMatch.Score"/>
        public double Score { get; }

        /// <inheritdoc cref="IHasSequenceCoverageFromFragments.FragmentCoveragePositionInPeptide"/>
        public HashSet<int>? FragmentCoveragePositionInPeptide { get; protected set; }

        /// <summary>
        /// N-terminal fragment positions (one-based). For peptides, these correspond to b-ion positions.
        /// For nucleic acids, these correspond to 5' fragment positions.
        /// Set via <see cref="SetFragmentPositions"/> before calling <see cref="GetSequenceCoverage"/>.
        /// </summary>
        protected List<int>? NTerminalFragmentPositions { get; set; }

        /// <summary>
        /// C-terminal fragment positions (one-based). For peptides, these correspond to y-ion positions.
        /// For nucleic acids, these correspond to 3' fragment positions.
        /// Set via <see cref="SetFragmentPositions"/> before calling <see cref="GetSequenceCoverage"/>.
        /// </summary>
        protected List<int>? CTerminalFragmentPositions { get; set; }

        private readonly List<IBioPolymerWithSetMods> _identifiedBioPolymers;

        /// <inheritdoc cref="ISpectralMatch.GetIdentifiedBioPolymersWithSetMods"/>
        public IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods()
        {
            return _identifiedBioPolymers;
        }

        /// <summary>
        /// Calculates which residues are covered by fragment ions and populates 
        /// <see cref="FragmentCoveragePositionInPeptide"/> with their one-based positions.
        /// 
        /// A residue is considered covered if:
        /// <list type="bullet">
        ///   <item><description>It has sequential N-terminal fragments on both sides</description></item>
        ///   <item><description>It has sequential C-terminal fragments on both sides</description></item>
        ///   <item><description>It has both N-terminal and C-terminal fragment evidence</description></item>
        ///   <item><description>It is the first or last residue with appropriate terminal fragment</description></item>
        /// </list>
        /// </summary>
        /// <remarks>
        /// Call <see cref="SetFragmentPositions"/> before this method to provide fragment data.
        /// If no fragment positions are set, <see cref="FragmentCoveragePositionInPeptide"/> remains null.
        /// This method can be called multiple times; each call recalculates coverage from the current fragment positions.
        /// </remarks>
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

            // Process N-terminal fragments (b-ions for peptides, 5' fragments for nucleic acids)
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

            // Process C-terminal fragments (y-ions for peptides, 3' fragments for nucleic acids)
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
        /// Compares to another <see cref="IHasSequenceCoverageFromFragments"/>.
        /// If the other object also implements <see cref="ISpectralMatch"/>, delegates to 
        /// <see cref="CompareTo(ISpectralMatch)"/>. Otherwise returns 0 (equal) or -1 if other is null.
        /// </summary>
        public int CompareTo(IHasSequenceCoverageFromFragments? other)
        {
            if (other is ISpectralMatch spectralMatch)
            {
                return CompareTo(spectralMatch);
            }
            return other is null ? -1 : 0;
        }

        /// <summary>
        /// Sets the fragment positions used by <see cref="GetSequenceCoverage"/> to calculate coverage.
        /// </summary>
        /// <param name="nTerminalPositions">One-based positions of N-terminal fragments (b-ions for peptides, 5' for nucleic acids).</param>
        /// <param name="cTerminalPositions">One-based positions of C-terminal fragments (y-ions for peptides, 3' for nucleic acids).</param>
        public void SetFragmentPositions(IEnumerable<int>? nTerminalPositions, IEnumerable<int>? cTerminalPositions)
        {
            NTerminalFragmentPositions = nTerminalPositions?.ToList();
            CTerminalFragmentPositions = cTerminalPositions?.ToList();
        }

        /// <summary>
        /// Adds a biopolymer to the list of identified biopolymers for this match.
        /// Null values are ignored.
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
        /// Null entries in the collection are filtered out.
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
        /// Orders by Score (descending), then FullFilePath (ascending), then OneBasedScanNumber (ascending).
        /// </summary>
        /// <returns>Negative if this should sort before other; positive if after; zero if equal.</returns>
        public int CompareTo(ISpectralMatch? other)
        {
            if (other is null) return -1;

            int scoreComparison = other.Score.CompareTo(Score); // Descending
            if (scoreComparison != 0) return scoreComparison;

            int fileComparison = string.Compare(FullFilePath, other.FullFilePath, StringComparison.Ordinal);
            if (fileComparison != 0) return fileComparison;

            return OneBasedScanNumber.CompareTo(other.OneBasedScanNumber);
        }

        /// <summary>
        /// Determines equality based on FullFilePath, OneBasedScanNumber, and FullSequence.
        /// Score is intentionally excluded from equality comparison.
        /// </summary>
        public bool Equals(SpectralMatchWithSequenceCoverage? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;

            return string.Equals(FullFilePath, other.FullFilePath, StringComparison.Ordinal)
                && OneBasedScanNumber == other.OneBasedScanNumber
                && string.Equals(FullSequence, other.FullSequence, StringComparison.Ordinal);
        }

        /// <inheritdoc/>
        public override bool Equals(object? obj)
        {
            if (obj is SpectralMatchWithSequenceCoverage sm) return Equals(sm);
            return false;
        }

        /// <inheritdoc/>
        public override int GetHashCode()
        {
            return HashCode.Combine(FullFilePath, OneBasedScanNumber, FullSequence);
        }

        /// <summary>
        /// Returns a human-readable string representation: "Scan {number}: {sequence} (Score: {score})".
        /// </summary>
        public override string ToString()
        {
            return $"Scan {OneBasedScanNumber}: {FullSequence} (Score: {Score:F2})";
        }

        /// <summary>
        /// Equality operator. Two matches are equal if they have the same FullFilePath, OneBasedScanNumber, and FullSequence.
        /// </summary>
        public static bool operator ==(SpectralMatchWithSequenceCoverage? left, SpectralMatchWithSequenceCoverage? right)
        {
            if (left is null) return right is null;
            return left.Equals(right);
        }

        /// <summary>
        /// Inequality operator.
        /// </summary>
        public static bool operator !=(SpectralMatchWithSequenceCoverage? left, SpectralMatchWithSequenceCoverage? right)
        {
            return !(left == right);
        }
    }
}