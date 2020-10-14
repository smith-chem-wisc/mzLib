namespace FlashLFQ
{
    public class SpectraFileInfo
    {
        public readonly string FullFilePathWithExtension;
        public readonly string FilenameWithoutExtension;
        public string SampleGroup;
        public readonly int Sample;
        public readonly int Fraction;
        public readonly int Replicate;

        public SpectraFileInfo(string fullFilePathWithExtension, string sampleGroup, int sample, int replicate, int fraction)
        {
            this.FullFilePathWithExtension = fullFilePathWithExtension;
            this.FilenameWithoutExtension = System.IO.Path.GetFileNameWithoutExtension(this.FullFilePathWithExtension);
            this.SampleGroup = sampleGroup;
            this.Sample = sample;
            this.Replicate = replicate;
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