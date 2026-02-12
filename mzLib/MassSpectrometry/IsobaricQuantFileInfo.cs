using System;
using System.Collections.Generic;
using System.IO;


namespace MassSpectrometry
{
    /// <summary>
    /// Metadata describing a single isobaric-quantification input file and its plex channel annotations.
    /// Instances encapsulate the file identity (path, plex, fraction, technical replicate) and an
    /// ordered, read-only set of per-channel annotations used for mapping reporter channels to samples.
    /// Designed as a small immutable value object for use as a dictionary key or in result tables.
    /// Note: file path identity is compared using a normalized (canonical) representation computed
    /// with <see cref="Path.GetFullPath(string)"/> when possible; if the provided path is null or empty
    /// it is normalized to an empty string.
    /// </summary>
    public class IsobaricQuantFileInfo : IEquatable<IsobaricQuantFileInfo>
    {
        private readonly string _canonicalFullFilePathWithExtension;

        /// <summary>
        /// Create a new instance describing a single isobaric quantification file.
        /// </summary>
        /// <param name="fullFilePathWithExtension">Absolute or relative path to the spectra/quant file including extension.
        /// If null or empty, it will be normalized to empty string.</param>
        /// <param name="plex">The plex identifier (for example "TMT10" or "TMTpro16"). Null is normalized to an empty string.</param>
        /// <param name="fraction">Fraction index (1-based) for fractionated experiments.</param>
        /// <param name="technicalReplicate">Technical replicate index (1-based) for repeated injections/runs.</param>
        /// <param name="annotations">Ordered list of channel annotations for the file's plex. May be null or empty.</param>
        public IsobaricQuantFileInfo(string fullFilePathWithExtension, string plex, int fraction, int technicalReplicate, IReadOnlyList<IsobaricQuantPlexAnnotation> annotations)
        {
            // Normalize visible properties
            FullFilePathWithExtension = fullFilePathWithExtension ?? string.Empty;
            Plex = plex ?? string.Empty;
            Fraction = fraction;
            TechnicalReplicate = technicalReplicate;
            Annotations = annotations ?? Array.Empty<IsobaricQuantPlexAnnotation>();

            // Compute a canonical path for identity comparisons.
            // Use Path.GetFullPath to resolve relative segments and normalize separators where possible.
            // If GetFullPath throws (invalid path characters etc.), fall back to the original (normalized) string.
            if (string.IsNullOrEmpty(FullFilePathWithExtension))
            {
                _canonicalFullFilePathWithExtension = string.Empty;
            }
            else
            {
                try
                {
                    _canonicalFullFilePathWithExtension = Path.GetFullPath(FullFilePathWithExtension);
                }
                catch
                {
                    // Keep the original string if canonicalization fails
                    _canonicalFullFilePathWithExtension = FullFilePathWithExtension;
                }
            }
        }

        /// <summary>
        /// Absolute or relative path to the data file (for example a .raw or .mzML), including extension.
        /// This property stores the original (normalized-null) input. For identity comparisons the
        /// canonicalized path computed at construction time is used.
        /// </summary>
        public string FullFilePathWithExtension { get; }

        /// <summary>
        /// The plex identifier for this file (for example "TMT10" or "TMTpro16").
        /// Returns an empty string when not provided.
        /// </summary>
        public string Plex { get; }

        /// <summary>
        /// Fraction index for the file. Uses 1-based indexing to match common experimental notation.
        /// </summary>
        public int Fraction { get; }             // 1-based

        /// <summary>
        /// Technical replicate index for the file. Uses 1-based indexing.
        /// </summary>
        public int TechnicalReplicate { get; }   // 1-based

        /// <summary>
        /// Ordered, read-only collection of per-channel annotations for this file's plex
        /// (one entry per reporter/tag). May be empty but never null.
        /// </summary>
        public IReadOnlyList<IsobaricQuantPlexAnnotation> Annotations { get; } // All tags for this file's plex

        /// <summary>
        /// Equality compares the canonical identity of the file info.
        /// Two <see cref="IsobaricQuantFileInfo"/> instances are considered equal when all of the following match:
        /// - canonicalized <see cref="FullFilePathWithExtension"/> (ordinal comparison),
        /// - <see cref="Plex"/> (ordinal comparison),
        /// - <see cref="Fraction"/>, 
        /// - <see cref="TechnicalReplicate"/>.
        /// Canonicalization uses <see cref="Path.GetFullPath(string)"/> where possible; if canonicalization
        /// failed during construction, the original path string is used.
        /// </summary>
        /// <param name="obj">Object to compare.</param>
        /// <returns>
        /// True when <paramref name="obj"/> is an <see cref="IsobaricQuantFileInfo"/> with the same
        /// canonical file path, plex, fraction and technical replicate; otherwise false.
        /// </returns>
        public override bool Equals(object obj)
        {
            return Equals(obj as IsobaricQuantFileInfo);
        }

        /// <summary>
        /// Type-safe equality comparison. See <see cref="Equals(object)"/>.
        /// </summary>
        /// <param name="other">Other instance to compare (may be null).</param>
        /// <returns>True when identity components match; otherwise false.</returns>
        public bool Equals(IsobaricQuantFileInfo other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;

            return StringComparer.Ordinal.Equals(_canonicalFullFilePathWithExtension ?? string.Empty, other._canonicalFullFilePathWithExtension ?? string.Empty)
                && StringComparer.Ordinal.Equals(Plex ?? string.Empty, other.Plex ?? string.Empty)
                && Fraction == other.Fraction
                && TechnicalReplicate == other.TechnicalReplicate;
        }

        /// <summary>
        /// Returns a stable hash code derived from the identity components used by <see cref="Equals(object)"/>.
        /// The hash is computed from:
        /// - canonicalized <see cref="FullFilePathWithExtension"/> (ordinal),
        /// - <see cref="Plex"/> (ordinal),
        /// - <see cref="Fraction"/>,
        /// - <see cref="TechnicalReplicate"/>.
        /// </summary>
        /// <returns>Hash code for the file identity.</returns>
        public override int GetHashCode()
        {
            int h1 = StringComparer.Ordinal.GetHashCode(_canonicalFullFilePathWithExtension ?? string.Empty);
            int h2 = StringComparer.Ordinal.GetHashCode(Plex ?? string.Empty);
            return HashCode.Combine(h1, h2, Fraction, TechnicalReplicate);
        }

        /// <summary>
        /// Equality operator overload to match value semantics.
        /// </summary>
        public static bool operator ==(IsobaricQuantFileInfo left, IsobaricQuantFileInfo right)
        {
            return Equals(left, right);
        }

        /// <summary>
        /// Inequality operator overload to match value semantics.
        /// </summary>
        public static bool operator !=(IsobaricQuantFileInfo left, IsobaricQuantFileInfo right)
        {
            return !Equals(left, right);
        }

        /// <summary>
        /// Returns the file name portion (without directory) of <see cref="FullFilePathWithExtension"/>.
        /// Suitable for display and logging.
        /// </summary>
        /// <returns>File name portion of the file path, or empty string when the path is null/empty.</returns>
        public override string ToString()
        {
            return Path.GetFileName(FullFilePathWithExtension ?? string.Empty);
        }
    }

    /// <summary>
    /// Annotation describing a single reporter/tag within a plex.
    /// Map a reporter channel to sample metadata: tag label, sample name, condition, and biological replicate.
    /// </summary>
    public class IsobaricQuantPlexAnnotation
    {
        /// <summary>
        /// Reporter tag identifier (for example "126", "127N", "128C" or reader-specific tag label).
        /// </summary>
        public string Tag { get; init; } = "";

        /// <summary>
        /// Human-readable sample name associated with this reporter channel.
        /// </summary>
        public string SampleName { get; init; } = "";

        /// <summary>
        /// Experimental condition or group for this reporter channel (for example "Control" or "Treatment").
        /// </summary>
        public string Condition { get; init; } = "";

        /// <summary>
        /// Biological replicate index for this channel (1-based). Zero or negative indicates 'not set'.
        /// </summary>
        public int BiologicalReplicate { get; init; }
    }
}
