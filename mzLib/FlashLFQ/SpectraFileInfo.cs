namespace FlashLFQ
{
    public class SpectraFileInfo
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
        {
            this.FullFilePathWithExtension = fullFilePathWithExtension;
            this.FilenameWithoutExtension = System.IO.Path.GetFileNameWithoutExtension(this.FullFilePathWithExtension);
            this.Condition = condition;
            this.BiologicalReplicate = biorep;
            this.TechnicalReplicate = techrep;
            this.Fraction = fraction;
        }

        // files are considered the same if the absolute file path is the same
        public override bool Equals(object obj)
        {
            return base.Equals(obj) && ((SpectraFileInfo)obj).FullFilePathWithExtension.Equals(this.FullFilePathWithExtension);
        }

        public override int GetHashCode()
        {
            return FullFilePathWithExtension.GetHashCode();
        }
    }
}