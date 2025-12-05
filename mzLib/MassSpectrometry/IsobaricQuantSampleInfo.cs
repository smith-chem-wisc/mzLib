using System;

namespace MassSpectrometry
{
    /// <summary>
    /// Concrete implementation of <see cref="IIsobaricQuantSampleInfo"/>.
    /// Represents sample information for isobaric (TMT/iTRAQ) quantification.
    /// </summary>
    public class IsobaricQuantSampleInfo : IIsobaricQuantSampleInfo
    {
        /// <inheritdoc />
        public string ChannelLabel { get; }

        /// <inheritdoc />
        public string Condition { get; }

        /// <inheritdoc />
        public int BiologicalReplicate { get; }

        /// <inheritdoc />
        public int TechnicalReplicate { get; }

        /// <inheritdoc />
        public int PlexId { get; }

        /// <inheritdoc />
        public int Fraction { get; }

        /// <inheritdoc />
        public double ReporterIonMz { get; }

        /// <inheritdoc />
        public string? SampleName { get; }

        /// <inheritdoc />
        public bool IsReferenceChannel { get; }

        /// <inheritdoc />
        public string DisplayName => SampleName ?? $"Plex{PlexId}_{ChannelLabel}";

        /// <inheritdoc />
        public string FullFilePathWithExtension { get; }

        /// <summary>
        /// Creates a new isobaric quantification sample info.
        /// </summary>
        /// <param name="channelLabel">The isobaric channel label (e.g., "126", "127N").</param>
        /// <param name="reporterIonMz">The reporter ion m/z value for this channel.</param>
        /// <param name="condition">The condition or experimental group.</param>
        /// <param name="biologicalReplicate">Biological replicate identifier.</param>
        /// <param name="technicalReplicate">Technical replicate identifier.</param>
        /// <param name="plexId">The plex identifier (default 1 for single-plex experiments).</param>
        /// <param name="fraction">Fraction identifier (default 0 for non-fractionated).</param>
        /// <param name="fullFilePathWithExtension">Full path to the source file.</param>
        /// <param name="sampleName">Optional sample name for reporting.</param>
        /// <param name="isReferenceChannel">True if this is a reference/normalization channel.</param>
        public IsobaricQuantSampleInfo(
            string channelLabel,
            double reporterIonMz,
            string condition,
            int biologicalReplicate = 1,
            int technicalReplicate = 1,
            int plexId = 1,
            int fraction = 0,
            string fullFilePathWithExtension = "",
            string? sampleName = null,
            bool isReferenceChannel = false)
        {
            ChannelLabel = channelLabel ?? throw new ArgumentNullException(nameof(channelLabel));
            ReporterIonMz = reporterIonMz;
            Condition = condition ?? string.Empty;
            BiologicalReplicate = biologicalReplicate;
            TechnicalReplicate = technicalReplicate;
            PlexId = plexId;
            Fraction = fraction;
            FullFilePathWithExtension = fullFilePathWithExtension ?? string.Empty;
            SampleName = sampleName;
            IsReferenceChannel = isReferenceChannel;
        }

        /// <inheritdoc />
        public bool Equals(IIsobaricQuantSampleInfo? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;
            return string.Equals(ChannelLabel, other.ChannelLabel, StringComparison.Ordinal)
                && PlexId == other.PlexId;
        }

        /// <summary>
        /// Override of Object.Equals to use the typed Equals.
        /// </summary>
        public override bool Equals(object? obj)
        {
            return Equals(obj as IIsobaricQuantSampleInfo);
        }

        /// <summary>
        /// Returns a hash code based on <see cref="ChannelLabel"/> and <see cref="PlexId"/>.
        /// </summary>
        public override int GetHashCode()
        {
            return HashCode.Combine(
                StringComparer.Ordinal.GetHashCode(ChannelLabel ?? string.Empty),
                PlexId);
        }

        /// <summary>
        /// Returns a string representation of this sample info.
        /// </summary>
        public override string ToString()
        {
            return DisplayName;
        }
    }
}