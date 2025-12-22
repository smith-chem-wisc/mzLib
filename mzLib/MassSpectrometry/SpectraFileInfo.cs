using MassSpectrometry.ExperimentalDesign;
using System;
using System.IO;

namespace MassSpectrometry
{
    /// <summary>
    /// Represents metadata for a spectra data file used in label-free quantification workflows.
    /// Uniquely identified by the combination of full file path, condition, biological replicate,
    /// technical replicate, and fraction. Used as a stable key across results, indexes, and normalization.
    /// Equality and hashing include all identity components to safely distinguish multiple entries
    /// that may share the same file path but differ by replicate/fraction.
    /// </summary>
    public class SpectraFileInfo : ISampleInfo, IEquatable<SpectraFileInfo>, IComparable<ISampleInfo>
    {
        /// <summary>
        /// The absolute path to the data file (e.g., a .raw or .mzML file) including the file extension.
        /// </summary>
        public string FullFilePathWithExtension { get; }

        /// <summary>
        /// The name of the data file without the directory path or file extension.
        /// Derived from <see cref="FullFilePathWithExtension"/> at construction time.
        /// </summary>
        public string FilenameWithoutExtension { get; }

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
        /// The fraction identifier for fractionated samples.
        /// Used when a single sample is separated into multiple fractions before analysis.
        /// Returns 0 if the sample is not fractionated.
        /// </summary>
        public int Fraction { get; }

        /// <summary>
        /// The technical replicate identifier for repeated instrument injections or runs.
        /// Distinguishes repeated measurements of the same biological sample.
        /// </summary>
        public int TechnicalReplicate { get; }

        /// <summary>
        /// Creates a new <see cref="SpectraFileInfo"/> instance with the specified parameters.
        /// </summary>
        /// <param name="fullFilePathWithExtension">The absolute path to the data file including extension.</param>
        /// <param name="condition">The sample group or condition label.</param>
        /// <param name="biorep">The biological replicate identifier.</param>
        /// <param name="techrep">The technical replicate identifier.</param>
        /// <param name="fraction">The fraction identifier (0 if not fractionated).</param>
        public SpectraFileInfo(string fullFilePathWithExtension, string condition, int biorep, int techrep, int fraction)
        {
            FullFilePathWithExtension = fullFilePathWithExtension;
            FilenameWithoutExtension = Path.GetFileNameWithoutExtension(FullFilePathWithExtension);
            Condition = condition;
            BiologicalReplicate = biorep;
            Fraction = fraction;
            TechnicalReplicate = techrep;
        }
        /// <summary>
        /// Determines whether two <see cref="SpectraFileInfo"/> instances are equal.
        /// Two samples are equal if they have the same <see cref="FullFilePathWithExtension"/> and <see cref="ChannelLabel"/>.
        /// </summary>
        /// <param name="left">The first instance to compare.</param>
        /// <param name="right">The second instance to compare.</param>
        /// <returns>True if both instances are equal or both are null; otherwise false.</returns>
        public static bool operator ==(SpectraFileInfo? left, SpectraFileInfo? right)
        {
            if (left is null) return right is null;
            return left.Equals(right);
        }

        /// <summary>
        /// Determines whether two <see cref="SpectraFileInfo"/> instances are not equal.
        /// Two samples are not equal if they differ in <see cref="FullFilePathWithExtension"/> or <see cref="ChannelLabel"/>.
        /// </summary>
        /// <param name="left">The first instance to compare.</param>
        /// <param name="right">The second instance to compare.</param>
        /// <returns>True if the instances are not equal; otherwise false.</returns>
        public static bool operator !=(SpectraFileInfo? left, SpectraFileInfo? right)
        {
            return !(left == right);
        }
        /// <summary>
        /// Determines whether the specified <see cref="SpectraFileInfo"/> is equal to this instance.
        /// Two instances are equal if they have the same file path, condition, biological replicate,
        /// fraction, and technical replicate.
        /// </summary>
        /// <param name="other">The other <see cref="SpectraFileInfo"/> to compare.</param>
        /// <returns>True if all identity components are equal; otherwise false.</returns>
        public bool Equals(SpectraFileInfo? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;

            return string.Equals(FullFilePathWithExtension, other.FullFilePathWithExtension, StringComparison.Ordinal)
                && string.Equals(Condition, other.Condition, StringComparison.Ordinal)
                && BiologicalReplicate == other.BiologicalReplicate
                && Fraction == other.Fraction
                && TechnicalReplicate == other.TechnicalReplicate;
        }

        /// <summary>
        /// Determines whether the specified <see cref="ISampleInfo"/> is equal to this instance.
        /// Returns true only if the other object is a <see cref="SpectraFileInfo"/> with matching identity components.
        /// </summary>
        /// <param name="other">The other <see cref="ISampleInfo"/> to compare.</param>
        /// <returns>True if equal; otherwise false.</returns>
        public bool Equals(ISampleInfo? other)
        {
            if (other is SpectraFileInfo specInfo)
            {
                return Equals(specInfo);
            }
            return false;
        }

        /// <summary>
        /// Determines whether the specified object is equal to this instance.
        /// </summary>
        /// <param name="obj">The object to compare.</param>
        /// <returns>True if the object is a <see cref="SpectraFileInfo"/> or <see cref="ISampleInfo"/> 
        /// with matching identity components; otherwise false.</returns>
        public override bool Equals(object? obj)
        {
            if (obj is null) return false;
            if (obj is SpectraFileInfo otherSpecInfo)
            {
                return Equals(otherSpecInfo);
            }
            if (obj is ISampleInfo otherSampleInfo)
            {
                return Equals(otherSampleInfo);
            }
            return false;
        }

        /// <summary>
        /// Returns a hash code based on all identity components: <see cref="FullFilePathWithExtension"/>,
        /// <see cref="Condition"/>, <see cref="BiologicalReplicate"/>, <see cref="Fraction"/>, 
        /// and <see cref="TechnicalReplicate"/>.
        /// </summary>
        /// <returns>A hash code for this instance.</returns>
        public override int GetHashCode()
        {
            return HashCode.Combine(
                StringComparer.Ordinal.GetHashCode(FullFilePathWithExtension ?? string.Empty),
                StringComparer.Ordinal.GetHashCode(Condition ?? string.Empty),
                BiologicalReplicate,
                Fraction,
                TechnicalReplicate);
        }

        /// <summary>
        /// Returns the file name (with extension) extracted from <see cref="FullFilePathWithExtension"/>.
        /// </summary>
        /// <returns>The file name portion of the full file path.</returns>
        public override string ToString()
        {
            return Path.GetFileName(FullFilePathWithExtension);
        }

        /// <summary>
        /// Compares this instance to another <see cref="ISampleInfo"/> for sorting purposes.
        /// Comparison order: <see cref="FullFilePathWithExtension"/>, <see cref="Condition"/>,
        /// <see cref="BiologicalReplicate"/>, <see cref="Fraction"/>, <see cref="TechnicalReplicate"/>.
        /// All comparisons use ascending order (A before B for strings, 1 before 2 for numbers).
        /// Non-null values sort before null.
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

            if (other is SpectraFileInfo) 
            {
                // Check for equality returns 0 if equal.
                if (Equals(other)) return 0;

                // Compare by Condition (A before B)
                int comparison = string.Compare(Condition ?? string.Empty, other.Condition ?? string.Empty, StringComparison.Ordinal);
                if (comparison != 0) return comparison;

                // Compare by BiologicalReplicate (1 before 2)
                comparison = BiologicalReplicate.CompareTo(other.BiologicalReplicate);
                if (comparison != 0) return comparison;

                // Compare by Fraction (1 before 2)
                comparison = Fraction.CompareTo(other.Fraction);
                if (comparison != 0) return comparison;

                // Compare by TechnicalReplicate (1 before 2)
                comparison = TechnicalReplicate.CompareTo(other.TechnicalReplicate);
				if (comparison != 0) return comparison;

                return string.Compare(FullFilePathWithExtension ?? string.Empty, other.FullFilePathWithExtension ?? string.Empty, StringComparison.Ordinal);
			}

            return 1; 
        }
    }
}