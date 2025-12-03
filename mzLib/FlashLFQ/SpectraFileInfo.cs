using System.IO;

namespace FlashLFQ
{
    // Temporary shim; remove after all call sites are updated.
    public class SpectraFileInfo : MassSpectrometry.SpectraFileInfo
    {
        /// <summary>
        /// The path to the data file (e.g., a .raw file) with the extension
        /// </summary>
        public string FullFilePathWithExtension { get; init; }
        /// <summary>
        /// The name of the data file without the extension
        /// </summary>
        public string FilenameWithoutExtension { get; init; }
        /// <summary>
        /// The condition of the sample (e.g., "Control" or "Treatment")
        /// </summary>
        public string Condition { get; set; }

        public int BiologicalReplicate { get; set; }
        public int TechnicalReplicate { get; set; }
        public int Fraction { get; set; }

        public SpectraFileInfo(string fullFilePathWithExtension, string condition, int biorep, int techrep, int fraction)
            : base(fullFilePathWithExtension, condition, biorep, techrep, fraction)
        {
        }
    }
}