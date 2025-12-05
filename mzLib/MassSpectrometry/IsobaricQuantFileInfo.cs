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
    /// </summary>
    public class IsobaricQuantFileInfo
    {
        /// <summary>
        /// Create a new instance describing a single isobaric quantification file.
        /// </summary>
        /// <param name="fullFilePathWithExtension">Absolute or relative path to the spectra/quant file including extension.</param>
        /// <param name="plex">The plex identifier (for example "TMT10" or "TMTpro16"). Null is normalized to an empty string.</param>
        /// <param name="fraction">Fraction index (1-based) for fractionated experiments.</param>
        /// <param name="technicalReplicate">Technical replicate index (1-based) for repeated injections/runs.</param>
        /// <param name="annotations">Ordered list of channel annotations for the file's plex. May be null or empty.</param>
        public IsobaricQuantFileInfo(string fullFilePathWithExtension, string plex, int fraction, int technicalReplicate, IReadOnlyList<IsobaricQuantPlexAnnotation> annotations)
        {
            FullFilePathWithExtension = fullFilePathWithExtension;
            Plex = plex ?? string.Empty;
            Fraction = fraction;
            TechnicalReplicate = technicalReplicate;
            Annotations = annotations ?? Array.Empty<IsobaricQuantPlexAnnotation>();
        }

        /// <summary>
        /// Absolute or relative path to the data file (for example a .raw or .mzML), including extension.
        /// This property forms the primary identity component used by <see cref="Equals(object)"/> and <see cref="GetHashCode"/>.
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
        /// Equality compares the canonical file path. Two <see cref="IsobaricQuantFileInfo"/>
        /// instances are considered equal when their <see cref="FullFilePathWithExtension"/>
        /// strings are equal (ordinal comparison).
        /// </summary>
        /// <param name="obj">Object to compare.</param>
        /// <returns>True when <paramref name="obj"/> is an <see cref="IsobaricQuantFileInfo"/> with the same file path; otherwise false.</returns>
        public override bool Equals(object obj)
        {
            return obj is IsobaricQuantFileInfo other && FullFilePathWithExtension.Equals(other.FullFilePathWithExtension);
        }

        /// <summary>
        /// Returns a stable hash code derived from <see cref="FullFilePathWithExtension"/> using ordinal string hashing.
        /// </summary>
        /// <returns>Hash code for the file identity.</returns>
        public override int GetHashCode()
        {
            return StringComparer.Ordinal.GetHashCode(FullFilePathWithExtension ?? string.Empty);
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
        public string Tag { get; set; } = "";

        /// <summary>
        /// Human-readable sample name associated with this reporter channel.
        /// </summary>
        public string SampleName { get; set; } = "";

        /// <summary>
        /// Experimental condition or group for this reporter channel (for example "Control" or "Treatment").
        /// </summary>
        public string Condition { get; set; } = "";

        /// <summary>
        /// Biological replicate index for this channel (1-based). Zero or negative indicates 'not set'.
        /// </summary>
        public int BiologicalReplicate { get; set; }
    }
}
