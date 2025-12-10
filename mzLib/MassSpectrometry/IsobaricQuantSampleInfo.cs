using System;

namespace MassSpectrometry
{
    /// <summary>
    /// Concrete implementation of <see cref="IIsobaricQuantSampleInfo"/>.
    /// Represents sample information for isobaric (TMT/iTRAQ) quantification.
    /// </summary>
    public class IsobaricQuantSampleInfo : IIsobaricQuantSampleInfo
    {
        /// <summary>
        /// Full path or identifier for the source file. May be empty for non-file-based samples.
        /// </summary>
        public string FullFilePathWithExtension { get; }

        /// <summary>
        /// Display name for the sample (used in headers and reports).
        /// </summary>
        public string SampleIdentifier { get; }

        /// <summary>
        /// The condition or experimental group this sample belongs to.
        /// </summary>
        public string Condition { get; }

        /// <summary>
        /// Biological replicate identifier within a condition.
        /// </summary>
        public int BiologicalReplicate { get; }

        /// <summary>
        /// Technical replicate identifier for repeated measurements.
        /// </summary>
        public int TechnicalReplicate { get; }

        /// <summary>
        /// Fraction identifier for fractionated workflows. Returns 0 if not applicable.
        /// </summary>
        public int Fraction { get; }

        /// <summary>
        /// The plex or multiplex set identifier. Distinguishes between different isobaric
        /// labeling experiments when multiple plexes are analyzed together (e.g., for bridge normalization).
        /// </summary>
        public int PlexId { get; }

        /// <summary>
        /// The isobaric channel label (e.g., "126", "127N", "127C", "128N" for TMT;
        /// "113", "114", "115" for iTRAQ). Serves as the primary identifier for the sample
        /// within a multiplexed experiment.
        /// </summary>
        public string ChannelLabel { get; }

        /// <summary>
        /// The reporter ion m/z value for this channel. Used to extract intensity values
        /// from MS2/MS3 spectra during quantification.
        /// </summary>
        public double ReporterIonMz { get; }

        /// <summary>
        /// True if this channel is used as a reference or normalization channel
        /// (e.g., pooled reference in TMT experiments).
        /// </summary>
        public bool IsReferenceChannel { get; }

        /// <summary>
        /// Creates a new isobaric quantification sample info.
        /// </summary>
        /// <param name="fullFilePathWithExtension">Full path to the source file.</param>
        /// <param name="sampleIdentifier">Display name for the sample.</param>
        /// <param name="condition">The condition or experimental group.</param>
        /// <param name="biologicalReplicate">Biological replicate identifier.</param>
        /// <param name="technicalReplicate">Technical replicate identifier.</param>
        /// <param name="fraction">Fraction identifier.</param>
        /// <param name="plexId">The plex identifier.</param>
        /// <param name="channelLabel">The isobaric channel label (e.g., "126", "127N").</param>
        /// <param name="reporterIonMz">The reporter ion m/z value for this channel.</param>
        /// <param name="isReferenceChannel">True if this is a reference/normalization channel.</param>
        public IsobaricQuantSampleInfo(
            string fullFilePathWithExtension,
            string sampleIdentifier,
            string condition,
            int biologicalReplicate,
            int technicalReplicate,
            int fraction,
            int plexId,
            string channelLabel,
            double reporterIonMz,
            bool isReferenceChannel)
        {
            FullFilePathWithExtension = fullFilePathWithExtension ?? string.Empty;
            SampleIdentifier = sampleIdentifier ?? string.Empty;
            Condition = condition ?? string.Empty;
            BiologicalReplicate = biologicalReplicate;
            TechnicalReplicate = technicalReplicate;
            Fraction = fraction;
            PlexId = plexId;
            ChannelLabel = channelLabel ?? throw new ArgumentNullException(nameof(channelLabel));
            ReporterIonMz = reporterIonMz;
            IsReferenceChannel = isReferenceChannel;
        }

        /// <summary>
        /// Determines whether the specified <see cref="IIsobaricQuantSampleInfo"/> is equal to this instance.
        /// Two samples are equal if they have the same channel label and plex ID.
        /// </summary>
        /// <param name="other">The other sample info to compare.</param>
        /// <returns>True if equal; otherwise false.</returns>
        public bool Equals(IIsobaricQuantSampleInfo? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;
            return string.Equals(ChannelLabel, other.ChannelLabel, StringComparison.Ordinal)
                && PlexId == other.PlexId;
        }

        /// <summary>
        /// Determines whether the specified object is equal to this instance.
        /// </summary>
        /// <param name="obj">The object to compare.</param>
        /// <returns>True if equal; otherwise false.</returns>
        public override bool Equals(object? obj)
        {
            return Equals(obj as IIsobaricQuantSampleInfo);
        }

        /// <summary>
        /// Returns a hash code based on <see cref="ChannelLabel"/> and <see cref="PlexId"/>.
        /// </summary>
        /// <returns>A hash code for this instance.</returns>
        public override int GetHashCode()
        {
            return HashCode.Combine(
                StringComparer.Ordinal.GetHashCode(ChannelLabel ?? string.Empty),
                PlexId);
        }

        /// <summary>
        /// Returns a string representation of this sample info.
        /// </summary>
        /// <returns>The sample identifier.</returns>
        public override string ToString()
        {
            return SampleIdentifier;
        }
    }
}