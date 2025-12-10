using System;

namespace MassSpectrometry
{
    /// <summary>
    /// Represents sample information for isobaric (TMT/iTRAQ) quantification.
    /// Extends <see cref="ISampleInfo"/> with isobaric-specific properties such as
    /// channel label, reporter ion m/z, and plex identifier.
    /// </summary>
    public class IsobaricQuantSampleInfo : ISampleInfo, IEquatable<IsobaricQuantSampleInfo>
    {
        /// <summary>
        /// Full path or identifier for the source file. May be empty for non-file-based samples.
        /// </summary>
        public string FullFilePathWithExtension { get; }

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
            string condition,
            int biologicalReplicate,
            int technicalReplicate,
            int fraction,
            int plexId,
            string channelLabel,
            double reporterIonMz,
            bool isReferenceChannel)
        {
            ChannelLabel = channelLabel;
            PlexId = plexId;
            FullFilePathWithExtension = fullFilePathWithExtension;
            Condition = condition;
            BiologicalReplicate = biologicalReplicate;
            TechnicalReplicate = technicalReplicate;
            Fraction = fraction;
            ReporterIonMz = reporterIonMz;
            IsReferenceChannel = isReferenceChannel;
        }

        /// <summary>
        /// Determines whether the specified <see cref="IsobaricQuantSampleInfo"/> is equal to this instance.
        /// Two samples are equal if they have the same <see cref="ChannelLabel"/> and <see cref="PlexId"/>.
        /// </summary>
        /// <param name="other">The other sample info to compare.</param>
        /// <returns>True if equal; otherwise false.</returns>
        public bool Equals(IsobaricQuantSampleInfo? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;
            return string.Equals(ChannelLabel, other.ChannelLabel, StringComparison.Ordinal)
                && PlexId == other.PlexId;
        }

        public bool Equals(ISampleInfo other)
        {
            if (other is IsobaricQuantSampleInfo otherIsobaric)
            {
                return Equals(otherIsobaric);
            }
            return false;
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (obj is IsobaricQuantSampleInfo otherIsobaric)
            {
                return Equals(otherIsobaric);
            }
            else if (obj is ISampleInfo otherSampleInfo)
            {
                return Equals(otherSampleInfo);
            }

            return false;
        }

        /// <summary>
        /// Returns a hash code based on <see cref="ChannelLabel"/> and <see cref="PlexId"/>.
        /// </summary>
        /// <returns>A hash code for this instance.</returns>
        public override int GetHashCode()
        {
            return HashCode.Combine(
                StringComparer.Ordinal.GetHashCode(ChannelLabel),
                PlexId);
        }

        /// <summary>
        /// Returns a string representation of this sample in the format "Plex{PlexId}_{ChannelLabel}".
        /// </summary>
        /// <returns>A display string for this sample.</returns>
        public override string ToString()
        {
            return $"Plex{PlexId}_{ChannelLabel}";
        }

        public int CompareTo(ISampleInfo other)
        {
            if (other is null) return -1;
            if (other is IsobaricQuantSampleInfo otherIsobaric)
            {
                int plexComparison = PlexId.CompareTo(otherIsobaric.PlexId);
                if (plexComparison != 0) return plexComparison;
                return string.Compare(ChannelLabel, otherIsobaric.ChannelLabel, StringComparison.Ordinal);
            }
                    return 1;
        }

    }
}