namespace Omics.SpectralMatch;

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
    /// If the given Spectral Match is a decoy
    /// </summary>
    bool IsDecoy { get; }

    /// <summary>
    /// The accession (unique identifier) of the identification
    /// </summary>
    public string Accession { get; }
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
    /// Numeric score for the match. Most implementations use a higher-is-better convention; callersBestMatchingBioPolymersWithSetMods
    /// should consult the implementer documentation if a different convention is used.
    /// The score is used for ranking, filtering and tie-breaking in comparison.
    /// </summary>
    double Score { get; }

    /// <summary>
    /// An array containing the intensities for the spectral match
    /// If no quantification wasn't performed, this will be null
    ///     if Multiple Reporter Ions are present, the order of intensities should match the order of the reporter ions
    ///     If LFQ was performed, this will be a single element array with the intensity value
    /// </summary>
    double[]? Intensities { get; }

    /// <summary>
    /// Returns the set (zero or more) of identified biopolymer objects (for example peptides or proteoforms
    /// with localized set modifications) associated with this spectral match. Implementations should return
    /// a stable enumeration (do not mutate the underlying collection while callers iterate) and document
    /// the ordering semantics (for example best-first).
    /// </summary>
    /// <returns>An enumerable of identified <see cref="IBioPolymerWithSetMods"/> instances; may be empty but not null.</returns>
    IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods();
}