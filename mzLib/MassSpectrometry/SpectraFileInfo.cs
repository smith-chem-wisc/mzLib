using System;
using System.IO;

namespace MassSpectrometry
{

    /// <summary>
    /// Uniquely identified by the combination of full file path, condition, biological replicate,
    /// technical replicate, and fraction. Used as a stable key across results, indexes, and normalization.
    /// Properties:
    /// - FullFilePathWithExtension: absolute path to the data file.
    /// - FilenameWithoutExtension: file name sans extension (derived).
    /// - Condition: sample group/condition label (e.g., Control/Treatment).
    /// - BiologicalReplicate: replicate identifier within a condition.
    /// - TechnicalReplicate: replicate identifier for repeated instrument injections/runs.
    /// - Fraction: fraction identifier for fractionated samples.
    /// Equality and hashing include all identity components to safely distinguish multiple entries
    /// that may share the same file path but differ by replicate/fraction.
    /// </summary>
    public class SpectraFileInfo : ISampleInfo
    {
        /// <summary>
        /// The path to the data file (e.g., a .raw file) with the extension.
        /// </summary>
        public string FullFilePathWithExtension { get; }

        /// <summary>
        /// The name of the data file without the extension.
        /// </summary>
        public string FilenameWithoutExtension { get; }

        /// <inheritdoc />
        public string Condition { get; }

        /// <inheritdoc />
        public int BiologicalReplicate { get; }

        /// <inheritdoc />
        public int TechnicalReplicate { get; }

        /// <inheritdoc />
        public int Fraction { get; }

        /// <inheritdoc />
        public string DisplayName => FilenameWithoutExtension;

        public SpectraFileInfo(string fullFilePathWithExtension, string condition, int biorep, int techrep, int fraction)
        {
            FullFilePathWithExtension = fullFilePathWithExtension ?? string.Empty;
            FilenameWithoutExtension = Path.GetFileNameWithoutExtension(FullFilePathWithExtension);
            Condition = condition ?? string.Empty;
            BiologicalReplicate = biorep;
            TechnicalReplicate = techrep;
            Fraction = fraction;
        }

        // Files are considered the same if all identity components are equal. Previously we used only filename. 
        // But, there is no rule about having multiple SpectraFileInfo objects with the same filename.
        public override bool Equals(object? obj)
        {
            if (obj is not SpectraFileInfo other) return false;

            return string.Equals(FullFilePathWithExtension, other.FullFilePathWithExtension, StringComparison.Ordinal)
                && string.Equals(Condition, other.Condition, StringComparison.Ordinal)
                && BiologicalReplicate == other.BiologicalReplicate
                && TechnicalReplicate == other.TechnicalReplicate
                && Fraction == other.Fraction;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(
                StringComparer.Ordinal.GetHashCode(FullFilePathWithExtension ?? string.Empty),
                StringComparer.Ordinal.GetHashCode(Condition ?? string.Empty),
                BiologicalReplicate,
                TechnicalReplicate,
                Fraction);
        }

        public override string ToString()
        {
            return Path.GetFileName(FullFilePathWithExtension);
        }
    }
}