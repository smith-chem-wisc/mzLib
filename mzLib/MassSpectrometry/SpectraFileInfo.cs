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
            FullFilePathWithExtension = fullFilePathWithExtension;
            FilenameWithoutExtension = Path.GetFileNameWithoutExtension(FullFilePathWithExtension);
            Condition = condition;
            BiologicalReplicate = biorep;
            TechnicalReplicate = techrep;
            Fraction = fraction;
        }

        // Files are considered the same if the absolute file path is the same.
        public override bool Equals(object? obj)
        {
            return obj is SpectraFileInfo other &&
                   string.Equals(other.FullFilePathWithExtension, FullFilePathWithExtension, System.StringComparison.Ordinal);
        }

        public override int GetHashCode()
        {
            return FullFilePathWithExtension?.GetHashCode() ?? 0;
        }

        public override string ToString()
        {
            return Path.GetFileName(FullFilePathWithExtension);
        }
    }
}