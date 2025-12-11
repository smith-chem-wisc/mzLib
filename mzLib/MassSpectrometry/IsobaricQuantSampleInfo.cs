using System;

namespace MassSpectrometry
{
    /// <summary>
    /// Represents sample information for isobaric (TMT/iTRAQ) quantification.
    /// Extends <see cref="ISampleInfo"/> with isobaric-specific properties such as
    /// channel label, reporter ion m/z, and plex identifier.
    /// 
    /// Unlike <see cref="SpectraFileInfo"/> which uses all identity components for equality,
    /// this class defines equality based solely on <see cref="PlexId"/> and <see cref="ChannelLabel"/>
    /// because these two properties uniquely identify a sample within an isobaric experiment.
    /// </summary>
    public class IsobaricQuantSampleInfo : ISampleInfo, IEquatable<IsobaricQuantSampleInfo>, IComparable<ISampleInfo>
    {
        /// <summary>
        /// The absolute path to the source data file including the file extension.
        /// May be empty for non-file-based samples or when the sample is defined
        /// independently of a specific spectra file.
        /// </summary>
        public string FullFilePathWithExtension { get; }

        /// <summary>
        /// The sample group or condition label (e.g., "Control", "Treatment").
        /// Used for grouping samples during normalization and statistical analysis.
        /// </summary>
        public string Condition { get; }

        /// <summary>
        /// The biological replicate identifier within a condition.
        /// Distinguishes independent biological samples within the same experimental condition.
        /// </summary>
        public int BiologicalReplicate { get; }

        /// <summary>
        /// The technical replicate identifier for repeated instrument injections or runs.
        /// Distinguishes repeated measurements of the same biological sample.
        /// </summary>
        public int TechnicalReplicate { get; }

        /// <summary>
        /// The fraction identifier for fractionated workflows.
        /// Used when a single sample is separated into multiple fractions before analysis.
        /// Returns 0 if the sample is not fractionated.
        /// </summary>
        public int Fraction { get; }

        /// <summary>
        /// The plex or multiplex set identifier. Distinguishes between different isobaric
        /// labeling experiments when multiple plexes are analyzed together (e.g., for bridge normalization).
        /// This is a primary identity component for equality comparisons.
        /// </summary>
        public int PlexId { get; }

        /// <summary>
        /// The isobaric channel label (e.g., "126", "127N", "127C", "128N" for TMT;
        /// "113", "114", "115" for iTRAQ). Serves as the primary identifier for the sample
        /// within a multiplexed experiment. This is a primary identity component for equality comparisons.
        /// </summary>
        public string ChannelLabel { get; }

        /// <summary>
        /// The reporter ion m/z value for this channel. Used to extract intensity values
        /// from MS2/MS3 spectra during quantification.
        /// </summary>
        public double ReporterIonMz { get; }

        /// <summary>
        /// True if this channel is used as a reference or normalization channel
        /// (e.g., pooled reference in TMT experiments). Reference channels are typically
        /// used for ratio-based normalization across plexes.
        /// </summary>
        public bool IsReferenceChannel { get; }

        /// <summary>
        /// A numeric identifier derived from <see cref="PlexId"/> and <see cref="ChannelLabel"/>.
        /// Computed at construction time for efficient lookups and comparisons.
        /// </summary>
        public int UniqueIdentifier { get; }

        /// <summary>
        /// Creates a new <see cref="IsobaricQuantSampleInfo"/> instance with the specified parameters.
        /// </summary>
        /// <param name="fullFilePathWithExtension">The absolute path to the source file including extension.</param>
        /// <param name="condition">The sample group or condition label.</param>
        /// <param name="biologicalReplicate">The biological replicate identifier.</param>
        /// <param name="technicalReplicate">The technical replicate identifier.</param>
        /// <param name="fraction">The fraction identifier (0 if not fractionated).</param>
        /// <param name="plexId">The plex or multiplex set identifier.</param>
        /// <param name="channelLabel">The isobaric channel label (e.g., "126", "127N").</param>
        /// <param name="reporterIonMz">The reporter ion m/z value for this channel.</param>
        /// <param name="isReferenceChannel">True if this is a reference or normalization channel.</param>
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
            UniqueIdentifier = HashCode.Combine(
                FullFilePathWithExtension,
                StringComparer.Ordinal.GetHashCode(ChannelLabel));
        }

        /// <summary>
        /// Determines whether two <see cref="IsobaricQuantSampleInfo"/> instances are equal.
        /// Two samples are equal if they have the same <see cref="PlexId"/> and <see cref="ChannelLabel"/>.
        /// </summary>
        /// <param name="left">The first instance to compare.</param>
        /// <param name="right">The second instance to compare.</param>
        /// <returns>True if both instances are equal or both are null; otherwise false.</returns>
        public static bool operator ==(IsobaricQuantSampleInfo? left, IsobaricQuantSampleInfo? right)
        {
            if (left is null) return right is null;
            return left.Equals(right);
        }

        /// <summary>
        /// Determines whether two <see cref="IsobaricQuantSampleInfo"/> instances are not equal.
        /// Two samples are not equal if they differ in <see cref="PlexId"/> or <see cref="ChannelLabel"/>.
        /// </summary>
        /// <param name="left">The first instance to compare.</param>
        /// <param name="right">The second instance to compare.</param>
        /// <returns>True if the instances are not equal; otherwise false.</returns>
        public static bool operator !=(IsobaricQuantSampleInfo? left, IsobaricQuantSampleInfo? right)
        {
            return !(left == right);
        }

        /// <summary>
        /// Determines whether the specified <see cref="IsobaricQuantSampleInfo"/> is equal to this instance.
        /// Two samples are considered equal if they have the same <see cref="PlexId"/> and <see cref="ChannelLabel"/>.
        /// 
        /// Note: Unlike <see cref="SpectraFileInfo"/> which compares all identity components,
        /// isobaric samples are uniquely identified by their plex and channel alone, as a single
        /// channel within a plex represents exactly one sample regardless of other metadata.
        /// </summary>
        /// <param name="other">The other <see cref="IsobaricQuantSampleInfo"/> to compare.</param>
        /// <returns>True if <see cref="PlexId"/> and <see cref="ChannelLabel"/> are equal; otherwise false.</returns>
        public bool Equals(IsobaricQuantSampleInfo? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;
            return string.Equals(FullFilePathWithExtension, other.FullFilePathWithExtension, StringComparison.Ordinal)
                && string.Equals(ChannelLabel, other.ChannelLabel, StringComparison.Ordinal);
        }

        /// <summary>
        /// Determines whether the specified <see cref="ISampleInfo"/> is equal to this instance.
        /// Returns true only if the other object is an <see cref="IsobaricQuantSampleInfo"/> 
        /// with matching <see cref="PlexId"/> and <see cref="ChannelLabel"/>.
        /// </summary>
        /// <param name="other">The other <see cref="ISampleInfo"/> to compare.</param>
        /// <returns>True if equal; otherwise false.</returns>
        public bool Equals(ISampleInfo? other)
        {
            if (other is IsobaricQuantSampleInfo otherIsobaric)
            {
                return Equals(otherIsobaric);
            }
            return false;
        }

        /// <summary>
        /// Determines whether the specified object is equal to this instance.
        /// Delegates to the typed <see cref="Equals(IsobaricQuantSampleInfo?)"/> method
        /// to ensure consistent equality semantics for both typed and untyped comparisons.
        /// </summary>
        /// <param name="obj">The object to compare.</param>
        /// <returns>True if the object is an <see cref="IsobaricQuantSampleInfo"/> or <see cref="ISampleInfo"/>
        /// with matching <see cref="PlexId"/> and <see cref="ChannelLabel"/>; otherwise false.</returns>
        public override bool Equals(object? obj)
        {
            if (obj is null) return false;
            if (obj is IsobaricQuantSampleInfo otherIsobaric)
            {
                return Equals(otherIsobaric);
            }
            if (obj is ISampleInfo otherSampleInfo)
            {
                return Equals(otherSampleInfo);
            }
            return false;
        }

        /// <summary>
        /// Returns a hash code based on <see cref="PlexId"/> and <see cref="ChannelLabel"/>,
        /// consistent with the equality comparison used in <see cref="Equals(IsobaricQuantSampleInfo?)"/>.
        /// </summary>
        /// <returns>A hash code for this instance.</returns>
        public override int GetHashCode()
        {
            return HashCode.Combine(
                StringComparer.Ordinal.GetHashCode(FullFilePathWithExtension ?? string.Empty),
                StringComparer.Ordinal.GetHashCode(ChannelLabel));
        }

        /// <summary>
        /// Returns a string representation of this sample in the format "Plex{PlexId}_{ChannelLabel}".
        /// </summary>
        /// <returns>A display string identifying this sample by plex and channel.</returns>
        public override string ToString()
        {
            return $"Plex{PlexId}_{ChannelLabel}";
        }

        /// <summary>
        /// Compares this instance to another <see cref="ISampleInfo"/> for sorting purposes.
        /// For <see cref="IsobaricQuantSampleInfo"/> instances, compares by <see cref="PlexId"/> first,
        /// then by <see cref="ChannelLabel"/>. All comparisons use ascending order 
        /// (1 before 2 for numbers, A before B for strings). Non-null values sort before null.
        /// 
        /// When comparing to a non-isobaric <see cref="ISampleInfo"/>, isobaric samples sort after
        /// other sample types.
        /// </summary>
        /// <param name="other">The other <see cref="ISampleInfo"/> to compare to.</param>
        /// <returns>
        /// A negative value if this instance precedes <paramref name="other"/>;
        /// zero if they are equal; a positive value if this instance follows <paramref name="other"/>.
        /// </returns>
        public int CompareTo(ISampleInfo? other)
        {
            // Non-null comes before null
            if (other is null) return -1;

            if (other is IsobaricQuantSampleInfo otherIsobaric)
            {
                // Compare by PlexId (1 before 2)
                int plexComparison = PlexId.CompareTo(otherIsobaric.PlexId);
                if (plexComparison != 0) return plexComparison;

                // Compare by ChannelLabel (A before B)
                return string.Compare(ChannelLabel, otherIsobaric.ChannelLabel, StringComparison.Ordinal);
            }

            // Isobaric samples sort after other sample types
            return 1;
        }
    }
}