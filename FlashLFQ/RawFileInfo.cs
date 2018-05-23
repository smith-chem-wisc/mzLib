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
        public MsDataFile dataFile;

        // condition/biorep/techrep/fraction info
        public readonly string condition;
        public readonly int biologicalReplicate;
        public readonly int fraction;
        public readonly int technicalReplicate;

        #endregion Public Fields

        #region Public Constructors

        public RawFileInfo(string fullFilePathWithExtension, string condition, int biorep, int techrep, int fraction)
        {
            this.fullFilePathWithExtension = fullFilePathWithExtension;
            this.filenameWithoutExtension = System.IO.Path.GetFileNameWithoutExtension(this.fullFilePathWithExtension);
            this.dataFile = null;
            clearAfterDone = true;
            this.condition = condition;
            this.biologicalReplicate = biorep;
            this.technicalReplicate = techrep;
            this.fraction = fraction;
        }

        public RawFileInfo(string fullFilePathWithExtension, string condition, int biorep, int techrep, int fraction, MsDataFile dataFile)
        {
            this.fullFilePathWithExtension = fullFilePathWithExtension;
            this.filenameWithoutExtension = System.IO.Path.GetFileNameWithoutExtension(this.fullFilePathWithExtension);
            this.dataFile = dataFile;
            clearAfterDone = false;
            this.condition = condition;
            this.biologicalReplicate = biorep;
            this.technicalReplicate = techrep;
            this.fraction = fraction;
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