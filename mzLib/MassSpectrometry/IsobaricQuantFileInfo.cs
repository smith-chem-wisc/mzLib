using System;
using System.Collections.Generic;
using System.IO;


namespace MassSpectrometry
{
    /// <summary>
    /// Describes an isobaric-quantification input file and the plex annotations for that file.
    /// Immutable for its core identity (path, plex, fraction, technical replicate and annotations).
    /// </summary>
    public class IsobaricQuantFileInfo
    {
        /// <summary>
        /// Create a new instance describing a single isobaric quantification file.
        /// </summary>
        /// <param name="fullFilePathWithExtension">Absolute or relative path to the spectra/quant file including extension.</param>
        /// <param name="plex">The plex identifier (e.g., "TMT10", "TMTpro16") for this file. Null will be normalized to empty string.</param>
        /// <param name="fraction">Fraction index (1-based) for fractionated experiments.</param>
        /// <param name="technicalReplicate">Technical replicate index (1-based) for repeated injections/runs.</param>
        /// <param name="annotations">Ordered list of channel annotations for the file's plex (may be empty).</param>
        public IsobaricQuantFileInfo(string fullFilePathWithExtension, string plex, int fraction, int technicalReplicate, IReadOnlyList<IsobaricQuantPlexAnnotation> annotations)
        {
            FullFilePathWithExtension = fullFilePathWithExtension;
            Plex = plex ?? string.Empty;
            Fraction = fraction;
            TechnicalReplicate = technicalReplicate;
            Annotations = annotations ?? Array.Empty<IsobaricQuantPlexAnnotation>();
        }

        /// <summary>
        /// The path to the data file (for example a .raw or .mzML) including the extension.
        /// Used as the primary identity for equality and hashing.
        /// </summary>
        public string FullFilePathWithExtension { get; }

        /// <summary>
        /// The plex identifier for this file (for example "TMT10" or "TMTpro16").
        /// Empty string when not provided.
        /// </summary>
        public string Plex { get; }

        /// <summary>
        /// Fraction index for the file. Uses 1-based indexing to match experimental notation.
        /// </summary>
        public int Fraction { get; }             // 1-based

        /// <summary>
        /// Technical replicate index for the file. Uses 1-based indexing.
        /// </summary>
        public int TechnicalReplicate { get; }   // 1-based

        /// <summary>
        /// All channel annotations for this file's plex (one entry per reporter/tag).
        /// Immutable reference to the provided read-only list (may be empty).
        /// </summary>
        public IReadOnlyList<IsobaricQuantPlexAnnotation> Annotations { get; } // All tags for this file's plex

        /// <summary>
        /// Returns true when the other object represents the same file path.
        /// Equality is based on <see cref="FullFilePathWithExtension"/>.
        /// </summary>
        /// <param name="obj">Other object to compare.</param>
        /// <returns>True when <paramref name="obj"/> is an <see cref="IsobaricQuantFileInfo"/> with the same file path.</returns>
        public override bool Equals(object obj)
        {
            if (base.Equals(obj))
            {
                return ((IsobaricQuantFileInfo)obj).FullFilePathWithExtension.Equals(FullFilePathWithExtension);
            }

            return false;
        }

        /// <summary>
        /// Returns a stable hash code based on <see cref="FullFilePathWithExtension"/>.
        /// </summary>
        /// <returns>Hash code.</returns>
        public override int GetHashCode()
        {
            return FullFilePathWithExtension.GetHashCode();
        }

        /// <summary>
        /// Returns the file name (without directory) for display and logging.
        /// </summary>
        /// <returns>File name portion of <see cref="FullFilePathWithExtension"/>.</returns>
        public override string ToString()
        {
            return Path.GetFileName(FullFilePathWithExtension);
        }
    }

    /// <summary>
    /// Annotation for a single isobaric reporter/tag within a plex.
    /// Used to map reporter channels to sample/condition/biological replicate metadata.
    /// </summary>
    public class IsobaricQuantPlexAnnotation
    {
        /// <summary>
        /// Reporter tag identifier (e.g., "126", "127N", "128C" or the channel label used by the reader).
        /// </summary>
        public string Tag { get; set; } = "";

        /// <summary>
        /// Human-readable sample name associated with this reporter channel.
        /// </summary>
        public string SampleName { get; set; } = "";

        /// <summary>
        /// Experimental condition or group for this reporter channel (e.g., "Control" or "Treatment").
        /// </summary>
        public string Condition { get; set; } = "";

        /// <summary>
        /// Biological replicate index for this channel (1-based). Zero or negative if not set.
        /// </summary>
        public int BiologicalReplicate { get; set; }
    }
}
