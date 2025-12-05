using System;
using System.Collections.Generic;

namespace Omics
{
    /// <summary>
    /// Represents a single spectrum-to-sequence identification (PSM) produced by a reader or search engine.
    /// Implementations provide identifying sequence strings, a numeric confidence score, the originating
    /// file identifier, and access to the biopolymer objects (peptides/proteoforms) that were assigned to
    /// the spectrum.  Implementers must provide consistent comparison and hashing semantics so matches can
    /// be ordered and used as dictionary/set keys.
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
        /// Typical implementations order by <see cref="Score"/> (higher => greater) and then apply
        /// deterministic tie-breakers (for example file/scan, sequence) to provide a stable ordering.
        /// </summary>
        /// <param name="other">Other spectral match to compare against (may be null).</param>
        /// <returns>Positive if this &gt; other, zero if equal, negative if this &lt; other.</returns>
        int CompareTo(ISpectralMatch? other);

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