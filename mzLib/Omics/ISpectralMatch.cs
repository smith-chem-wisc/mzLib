namespace Omics
{
    /// <summary>
    /// Represents a single spectrum-to-sequence identification (PSM) produced by a reader or search engine.
    /// Implementations provide identifying sequence strings, a numeric confidence score, the originating
    /// file identifier, and access to the biopolymer objects (peptides/proteoforms) that were assigned to
    /// the spectrum. Implementations must provide a consistent comparison implementation so matches can
    /// be ordered when needed.
    /// </summary>
    public interface ISpectralMatch : IComparable<ISpectralMatch>
    {
        /// <summary>
        /// The (absolute or relative) file path or file identifier for the spectra file that produced this match.
        /// Implementations may normalize this value (for example to an empty string) but should document whether
        /// comparisons treat path differences (case or separators) specially.
        /// </summary>
        string FullFilePath { get; }

        /// <summary>
        /// The full modified sequence string in the reader/search format (includes modification annotations).
        /// This is the canonical representation produced by the parser/search and should be used when
        /// presenting or comparing modified sequences.
        /// </summary>
        string FullSequence { get; }

        /// <summary>
        /// The base (unmodified) sequence. Use this for sequence-level grouping and comparisons when
        /// modification details are not relevant.
        /// </summary>
        string BaseSequence { get; }

        /// <summary>
        /// The scan number associated with this identification.
        /// Convention is one-based indexing; implementers should document the indexing they use.
        /// </summary>
        int OneBasedScanNumber { get; }

        /// <summary>
        /// Numeric score for the match. Most implementations use a higher-is-better convention; callers
        /// should consult the implementer documentation if a different convention is used.
        /// The score is used for ranking, filtering and tie-breaking in comparison.
        /// </summary>
        double Score { get; }

        /// <summary>
        /// Positions in the biopolymer sequence (one-based) that are covered by fragment ions.
        /// Populated by <see cref="GetSequenceCoverage"/> when fragment coverage data is available.
        /// May be null if coverage has not been calculated or is not available.
        /// </summary>
        HashSet<int>? FragmentCoveragePositionInPeptide { get; }

        /// <summary>
        /// Returns the set (zero or more) of identified biopolymer objects (for example peptides or proteoforms
        /// with localized set modifications) associated with this spectral match. Implementations should return
        /// a stable enumeration (do not mutate the underlying collection while callers iterate) and document
        /// the ordering semantics (for example best-first).
        /// </summary>
        /// <returns>An enumerable of identified <see cref="IBioPolymerWithSetMods"/> instances; may be empty but not null.</returns>
        IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods();

        /// <summary>
        /// Calculates sequence coverage from fragment ions for this spectral match.
        /// Populates <see cref="FragmentCoveragePositionInPeptide"/> with one-based positions
        /// of residues that are covered by matched fragment ions.
        /// Works for any biopolymer type (proteins, nucleic acids, etc.).
        /// Implementations without fragment ion data may leave the property null or empty.
        /// </summary>
        void GetSequenceCoverage();
    }
}