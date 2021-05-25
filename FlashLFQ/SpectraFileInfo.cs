namespace FlashLFQ
{
    public class SpectraFileInfo
    {
        public readonly string FullFilePathWithExtension;
        public readonly string FilenameWithoutExtension;
        public string Condition;
        public readonly int BiologicalReplicate;
        public readonly int Fraction;
        public readonly int TechnicalReplicate;

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