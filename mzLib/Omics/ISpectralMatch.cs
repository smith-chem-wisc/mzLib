using System;
using System.Collections.Generic;

namespace Omics
{
    /// <summary>
    /// Representation of a spectral match (identification/hypothesis) produced by a reader/parser.
    /// Implementations represent a single PSM / spectral identification and provide:
    /// - identifying sequences (full modified and base),
    /// - a scored ranking (Score),
    /// - and a set of candidate hypotheses that explain the spectrum.
    /// Implementers should provide consistent comparison and hashing semantics so matches can be
    /// ordered and used as dictionary/set keys.
    /// </summary>
    public interface ISpectralMatch : IComparable<ISpectralMatch>
    {
        /// <summary>
        /// Absolute (or original) file path for the spectra file that produced this match.
        /// Implementations may normalize this value (for example to an empty string) but should
        /// document how equality / hashing treats file path differences.
        /// </summary>
        string FullFilePath { get; init; }

        /// <summary>
        /// Full modified sequence in the source reader's format (includes modification annotations).
        /// This should be the canonical, reader-provided representation of the identified sequence.
        /// </summary>
        string FullSequence { get; init; }

        /// <summary>
        /// Base sequence (no modifications). Use this for sequence-based grouping and comparisons
        /// when modification detail is not relevant.
        /// </summary>
        string BaseSequence { get; init; }

        /// <summary>
        /// Collection of best-matching hypotheses (candidate biopolymers with set modifications)
        /// produced by the search/assignment step. Implementations should populate this list in
        /// descending confidence order or document the ordering semantics.
        /// </summary>
        List<ISpectralMatchHypothesis> BestMatchingBioPolymersWithSetMods { get; init; }

        /// <summary>
        /// Numeric score for this spectral match. Semantics (higher-is-better or lower-is-better)
        /// should be documented by the implementer; most callers expect higher scores to indicate
        /// better matches.
        /// </summary>
        double Score { get; init; }

        /// <summary>
        /// Compare this match to another. Implementations should return:
        /// - a positive value if this match is greater than <paramref name="other"/>,
        /// - zero if they are equal,
        /// - a negative value if this match is less than <paramref name="other"/>.
        /// Typical implementations order by Score (higher => greater) and break ties deterministically
        /// (for example by file/scan or sequence) to provide a stable ordering.
        /// This method corresponds to IComparable&lt;ISpectralMatch&gt;.CompareTo.
        /// </summary>
        int CompareTo();

        /// <summary>
        /// Return a hash code suitable for use in dictionaries and hash sets. Implementations should
        /// compute a stable hash based on the identity of the match (for example file + scan + sequence
        /// + charge) so that equal matches produce equal hash codes.
        /// </summary>
        int GetHashCode();
    }
}