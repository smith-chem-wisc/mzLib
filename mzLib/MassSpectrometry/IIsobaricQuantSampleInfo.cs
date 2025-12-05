using System;

namespace MassSpectrometry
{
    /// <summary>
    /// Represents sample information for isobaric (TMT/iTRAQ) quantification.
    /// Each instance corresponds to a single sample within an isobaric labeling experiment,
    /// identified by the combination of channel, condition, and replicate information.
    /// Extends <see cref="ISampleInfo"/> with isobaric-specific properties.
    /// </summary>
    public interface IIsobaricQuantSampleInfo : ISampleInfo, IEquatable<IIsobaricQuantSampleInfo>
    {
        /// <summary>
        /// The isobaric channel label (e.g., "126", "127N", "127C", "128N" for TMT;
        /// "113", "114", "115" for iTRAQ). Serves as the primary identifier for the sample
        /// within a multiplexed experiment.
        /// </summary>
        string ChannelLabel { get; }

        /// <summary>
        /// The plex or multiplex set identifier. Distinguishes between different isobaric
        /// labeling experiments when multiple plexes are analyzed together (e.g., for bridge normalization).
        /// </summary>
        int PlexId { get; }

        /// <summary>
        /// The reporter ion m/z value for this channel. Used to extract intensity values
        /// from MS2/MS3 spectra during quantification.
        /// </summary>
        double ReporterIonMz { get; }

        /// <summary>
        /// Optional sample name or description for reporting purposes.
        /// </summary>
        string SampleName { get; }

        /// <summary>
        /// True if this channel is used as a reference or normalization channel
        /// (e.g., pooled reference in TMT experiments).
        /// </summary>
        bool IsReferenceChannel { get; }

        /// <summary>
        /// Default equality implementation based on channel label and plex ID.
        /// Two samples are equal if they represent the same channel in the same plex.
        /// </summary>
        bool IEquatable<IIsobaricQuantSampleInfo>.Equals(IIsobaricQuantSampleInfo? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;
            return string.Equals(ChannelLabel, other.ChannelLabel, StringComparison.Ordinal)
                && PlexId == other.PlexId;
        }
    }
}