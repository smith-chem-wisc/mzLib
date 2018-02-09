using MassSpectrometry;

namespace FlashLFQ
{
    public class RawFileInfo
    {
        #region Public Fields

        public readonly string fullFilePathWithExtension;
        public readonly string filenameWithoutExtension;
        public readonly bool clearAfterDone;
        public string analysisSummary;
        public IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> dataFile;

        #endregion Public Fields

        #region Public Constructors

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

        #endregion Public Constructors

        #region Public Methods

        // files are considered the same if the absolute file path is the same
        public override bool Equals(object obj)
        {
            return base.Equals(obj) && ((RawFileInfo)obj).fullFilePathWithExtension.Equals(this.fullFilePathWithExtension);
        }

        public override int GetHashCode()
        {
            return fullFilePathWithExtension.GetHashCode();
        }

        #endregion Public Methods
    }
}