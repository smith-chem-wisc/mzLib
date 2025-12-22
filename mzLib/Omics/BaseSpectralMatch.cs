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
    public class BaseSpectralMatch : ISpectralMatch, IHasSequenceCoverageFromFragments, IEquatable<ISpectralMatch>
    {
        private readonly List<IBioPolymerWithSetMods> _identifiedBioPolymers;

        /// <summary>
        /// Creates a new spectral match with the specified properties.
        /// </summary>
        /// <param name="fullFilePath">The file path or identifier for the source spectra file. Null values are converted to empty string.</param>
        /// <param name="oneBasedScanNumber">The one-based scan number for this identification.</param>
        /// <param name="score">The numeric score for this match (higher is better by convention).</param>
        /// <param name="fullSequence">The full modified sequence string with modification annotations. Null values are converted to empty string.</param>
        /// <param name="baseSequence">The unmodified base sequence. Null values are converted to empty string.</param>
        /// <param name="identifiedBioPolymers">Optional biopolymers identified for this match. A defensive copy is made.</param>
        public BaseSpectralMatch(
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

        #region ISpectralMatch Properties

        /// <summary>
        /// The file path or identifier for the source spectra file.
        /// Used to associate this match with its originating data file.
        /// </summary>
        /// <inheritdoc cref="ISpectralMatch.FullFilePath"/>
        public string FullFilePath { get; }

        /// <summary>
        /// The full modified sequence string with modification annotations.
        /// For peptides, this includes bracket notation for modifications (e.g., "PEP[Phospho]TIDE").
        /// For nucleic acids, this includes modification annotations in the appropriate format.
        /// </summary>
        /// <inheritdoc cref="ISpectralMatch.FullSequence"/>
        public string FullSequence { get; }

        /// <summary>
        /// The unmodified base sequence containing only the residue letters.
        /// For peptides, this is the amino acid sequence without any modification annotations.
        /// For nucleic acids, this is the nucleotide sequence.
        /// </summary>
        /// <inheritdoc cref="ISpectralMatch.BaseSequence"/>
        public string BaseSequence { get; }

        /// <summary>
        /// The one-based scan number for this identification.
        /// This corresponds to the scan index in the source spectra file (first scan = 1).
        /// </summary>
        /// <inheritdoc cref="ISpectralMatch.OneBasedScanNumber"/>
        public int OneBasedScanNumber { get; }

        /// <summary>
        /// The numeric score for this match. Higher values indicate better matches by convention.
        /// The specific scoring algorithm and scale depend on the search engine used.
        /// </summary>
        /// <inheritdoc cref="ISpectralMatch.Score"/>
        public double Score { get; }

        public double[]? QuantValues { get; set; }

        #endregion

        #region IHasSequenceCoverageFromFragments Properties

        /// <summary>
        /// One-based positions of residues covered by matched fragment ions.
        /// Null until <see cref="GetSequenceCoverage(IEnumerable{int}, IEnumerable{int})"/> is called.
        /// A residue is considered covered if it has fragment evidence from sequential ions
        /// or from both N-terminal and C-terminal directions.
        /// </summary>
        /// <inheritdoc cref="IHasSequenceCoverageFromFragments.FragmentCoveragePositionInPeptide"/>
        public HashSet<int>? FragmentCoveragePositionInPeptide { get; protected set; }

        #endregion

        #region Methods

        /// <summary>
        /// Returns the biopolymers (peptides, oligonucleotides, etc.) identified for this spectral match.
        /// For ambiguous matches, this may return multiple candidates.
        /// The returned enumerable is a snapshot of the current identifications.
        /// </summary>
        /// <returns>An enumerable of identified biopolymers with their modifications.</returns>
        /// <inheritdoc cref="ISpectralMatch.GetIdentifiedBioPolymersWithSetMods"/>
        public IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods()
        {
            return _identifiedBioPolymers;
        }

        /// <summary>
        /// Calculates which residues are covered by fragment ions and populates 
        /// <see cref="FragmentCoveragePositionInPeptide"/> with their one-based positions.
        /// This parameterless overload does nothing - use <see cref="GetSequenceCoverage(IEnumerable{int}, IEnumerable{int})"/>
        /// or override in derived classes that have access to fragment ion data.
        /// </summary>
        /// <remarks>
        /// This method exists to satisfy the <see cref="IHasSequenceCoverageFromFragments"/> interface.
        /// Derived classes that store matched fragment ions should override this method to extract
        /// fragment positions and delegate to <see cref="GetSequenceCoverage(IEnumerable{int}, IEnumerable{int})"/>.
        /// </remarks>
        public virtual void GetSequenceCoverage()
        {
            // Base implementation does nothing - derived classes or callers should use
            // the overload that accepts fragment positions, or override this method
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
        /// <param name="nTerminalPositions">One-based positions of N-terminal fragments (b-ions for peptides, 5' for nucleic acids). Null is treated as empty.</param>
        /// <param name="cTerminalPositions">One-based positions of C-terminal fragments (y-ions for peptides, 3' for nucleic acids). Null is treated as empty.</param>
        /// <remarks>
        /// This method can be called multiple times; each call recalculates coverage from the provided fragment positions.
        /// If <see cref="BaseSequence"/> is null or empty, this method returns without modifying coverage.
        /// If both position lists are empty, this method returns without modifying coverage.
        /// </remarks>
        public void GetSequenceCoverage(IEnumerable<int>? nTerminalPositions, IEnumerable<int>? cTerminalPositions)
        {
            if (string.IsNullOrEmpty(BaseSequence))
            {
                return;
            }

            var nTermPositions = nTerminalPositions?.ToList() ?? new List<int>();
            var cTermPositions = cTerminalPositions?.ToList() ?? new List<int>();

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
        /// Adds a biopolymer to the list of identified biopolymers for this match.
        /// Use this method when additional biopolymer candidates are identified after construction.
        /// </summary>
        /// <param name="bioPolymer">The biopolymer to add. Null values are ignored.</param>
        public void AddIdentifiedBioPolymer(IBioPolymerWithSetMods bioPolymer)
        {
            if (bioPolymer != null)
            {
                _identifiedBioPolymers.Add(bioPolymer);
            }
        }

        /// <summary>
        /// Adds multiple biopolymers to the list of identified biopolymers for this match.
        /// Use this method when additional biopolymer candidates are identified after construction.
        /// </summary>
        /// <param name="bioPolymers">The biopolymers to add. Null entries in the collection are filtered out.</param>
        public void AddIdentifiedBioPolymers(IEnumerable<IBioPolymerWithSetMods> bioPolymers)
        {
            if (bioPolymers != null)
            {
                _identifiedBioPolymers.AddRange(bioPolymers.Where(b => b != null));
            }
        }

        #endregion

        #region IComparable Implementation

        /// <summary>
        /// Compares this spectral match to another for ordering purposes.
        /// Orders by Score (descending), then FullFilePath (ascending), then OneBasedScanNumber (ascending).
        /// This ordering ensures higher-scoring matches appear first, with ties broken deterministically.
        /// </summary>
        /// <param name="other">The other spectral match to compare to.</param>
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

        public bool Equals(ISpectralMatch? other)
        {
            if (other is BaseSpectralMatch sm) return Equals(sm);
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