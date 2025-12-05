

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
        int ScanNumber { get; }

        /// <summary>
        /// Numeric score for the match. Most implementations use a higher-is-better convention; callers
        /// should consult the implementer documentation if a different convention is used.
        /// The score is used for ranking, filtering and tie-breaking in comparison.
        /// </summary>
        double Score { get; }

        /// <summary>
        /// Compare this match to <paramref name="other"/> for ordering.
        /// Implementations should:
        /// - return a positive value when this match is greater than <paramref name="other"/>,
        /// - return zero when they are considered equal,
        /// - return a negative value when this match is less than <paramref name="other"/>.
        /// The default implementation orders by:
        /// 1) <see cref="Score"/> descending (higher scores are considered greater),
        /// 2) <see cref="FullFilePath"/> ascending (ordinal),
        /// 3) <see cref="FullSequence"/> ascending (ordinal),
        /// 4) <see cref="BaseSequence"/> ascending (ordinal),
        /// 5) <see cref="ScanNumber"/> ascending.
        /// Implementers may override to provide different tie-breaking semantics.
        /// </summary>
        /// <param name="other">Other spectral match to compare against (may be null).</param>
        /// <returns>Positive if this &gt; other, zero if equal, negative if this &lt; other.</returns>
        int CompareTo(ISpectralMatch? other)
        {
            if (other is null) return 1;

            // Primary: Score (higher is better) -> descending order
            int scoreCmp = Score.CompareTo(other.Score);
            if (scoreCmp != 0) return scoreCmp;

            // Tie-breakers: ascending order (ordinal)
            int fileCmp = string.Compare(FullFilePath ?? string.Empty, other.FullFilePath ?? string.Empty, StringComparison.Ordinal);
            if (fileCmp != 0) return fileCmp;

            int fullSeqCmp = string.Compare(FullSequence ?? string.Empty, other.FullSequence ?? string.Empty, StringComparison.Ordinal);
            if (fullSeqCmp != 0) return fullSeqCmp;

            int baseSeqCmp = string.Compare(BaseSequence ?? string.Empty, other.BaseSequence ?? string.Empty, StringComparison.Ordinal);
            if (baseSeqCmp != 0) return baseSeqCmp;

            // Final tie-breaker: scan number ascending
            int scanCmp = ScanNumber.CompareTo(other.ScanNumber);
            if (scanCmp != 0) return scanCmp;

            return 0;
        }

        /// <summary>
        /// Returns the set (zero or more) of identified biopolymer objects (for example peptides or proteoforms
        /// with localized set modifications) associated with this spectral match. Implementations should return
        /// a stable enumeration (do not mutate the underlying collection while callers iterate) and document
        /// the ordering semantics (for example best-first).
        /// </summary>
        /// <returns>An enumerable of identified <see cref="IBioPolymerWithSetMods"/> instances; may be empty but not null.</returns>
        IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods();
    }
}