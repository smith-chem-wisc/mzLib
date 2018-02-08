using System.Collections.Generic;
using System.Text;

namespace FlashLFQ
{
    public class ProteinGroup
    {
        #region Public Fields

        public static List<RawFileInfo> rawFiles;
        public readonly string ProteinGroupName;
        public Dictionary<RawFileInfo, double> intensities;
        
        #endregion Public Fields

        #region Public Properties

        public static string TabSeparatedHeader
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append("Protein Groups" + "\t");
                foreach (var rawfile in rawFiles)
                    sb.Append("Intensity_" + rawfile.filenameWithoutExtension + "\t");
                return sb.ToString();
            }
        }

        #endregion Public Properties

        #region Public Constructors

        public ProteinGroup(string proteinGroupName)
        {
            this.ProteinGroupName = proteinGroupName;
            this.intensities = new Dictionary<RawFileInfo, double>();

            foreach (var file in rawFiles)
                intensities.Add(file, 0);
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            StringBuilder str = new StringBuilder();
            str.Append(ProteinGroupName + "\t");

            foreach (var file in rawFiles)
                str.Append(intensities[file] + "\t");

            return str.ToString();
        }

        #endregion Public Methods
    }
}