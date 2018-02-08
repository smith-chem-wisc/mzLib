using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Text;

namespace FlashLFQ
{
    public class RawFileInfo
    {
        public readonly string fullFilePathWithExtension;
        public readonly string filenameWithoutExtension;
        public string analysisSummary;
        public readonly bool clearAfterDone;
        public IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> dataFile;

        public RawFileInfo(string fullFilePathWithExtension)
        {
            this.fullFilePathWithExtension = fullFilePathWithExtension;
            this.filenameWithoutExtension = System.IO.Path.GetFileNameWithoutExtension(this.fullFilePathWithExtension);
            this.dataFile = null;
            clearAfterDone = true;
        }

        public RawFileInfo(string fullFilePathWithExtension, IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> dataFile)
        {
            this.fullFilePathWithExtension = fullFilePathWithExtension;
            this.filenameWithoutExtension = System.IO.Path.GetFileNameWithoutExtension(this.fullFilePathWithExtension);
            this.dataFile = dataFile;
            clearAfterDone = false;
        }

        // files are considered the same if the absolute file path is the same
        public override bool Equals(object obj)
        {
            return base.Equals(obj) && ((RawFileInfo) obj).fullFilePathWithExtension.Equals(this.fullFilePathWithExtension);
        }

        public override int GetHashCode()
        {
            return fullFilePathWithExtension.GetHashCode();
        }
    }
}
