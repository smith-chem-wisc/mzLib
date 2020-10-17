using System.Collections.Generic;

namespace FlashLFQ
{
    public class SpectraFileInfo
    {
        public readonly string FullFilePathWithExtension;
        public readonly string FilenameWithoutExtension;
        public readonly List<Sample> Samples;

        public SpectraFileInfo(string fullFilePathWithExtension)
        {
            this.FullFilePathWithExtension = fullFilePathWithExtension;
            this.FilenameWithoutExtension = System.IO.Path.GetFileNameWithoutExtension(this.FullFilePathWithExtension);
            this.Samples = new List<Sample>();
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