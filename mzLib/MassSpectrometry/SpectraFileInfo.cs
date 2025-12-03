using System;
using System.IO;

namespace MassSpectrometry
{
    // Immutable value object representing a spectra file and its sample metadata.
    public sealed class SpectraFileInfo
    {
        /// <summary>
        /// The path to the data file (e.g., a .raw file) with the extension.
        /// </summary>
        public string FullFilePathWithExtension { get; }

        /// <summary>
        /// The name of the data file without the extension.
        /// </summary>
        public string FilenameWithoutExtension { get; }

        /// <summary>
        /// The condition of the sample (e.g., "Control" or "Treatment").
        /// </summary>
        public string Condition { get; }

        public int BiologicalReplicate { get; }
        public int TechnicalReplicate { get; }
        public int Fraction { get; }

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