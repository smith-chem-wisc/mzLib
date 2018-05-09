using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FlashLFQ
{
    public enum DetectionType { MSMS, MBR, NotDetected, MSMSAmbiguousPeakfinding, MSMSIdentifiedButNotQuantified }

    public class Peptide
    {
        #region Public Fields

        public static List<RawFileInfo> rawFiles;
        public readonly string Sequence;
        public Dictionary<RawFileInfo, double> intensities;
        public Dictionary<RawFileInfo, DetectionType> detectionTypes;
        public HashSet<ProteinGroup> proteinGroups;

        #endregion Public Fields

        #region Public Constructors

        public Peptide(string sequence)
        {
            this.Sequence = sequence;
            intensities = new Dictionary<RawFileInfo, double>();
            detectionTypes = new Dictionary<RawFileInfo, DetectionType>();
            proteinGroups = new HashSet<ProteinGroup>();

            foreach (var file in rawFiles)
                intensities.Add(file, 0);
            foreach (var file in rawFiles)
                detectionTypes.Add(file, DetectionType.NotDetected);
        }

        #endregion Public Constructors

        #region Public Properties

        public static string TabSeparatedHeader
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append("Sequence" + "\t");
                sb.Append("Protein Groups" + "\t");
                sb.Append("Gene Names" + "\t");
                sb.Append("Organism" + "\t");
                foreach (var rawfile in rawFiles)
                    sb.Append("Intensity_" + rawfile.filenameWithoutExtension + "\t");
                foreach (var rawfile in rawFiles)
                    sb.Append("Detection Type_" + rawfile.filenameWithoutExtension + "\t");
                return sb.ToString();
            }
        }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            StringBuilder str = new StringBuilder();
            str.Append(Sequence + "\t");
            str.Append(string.Join(";", proteinGroups.Select(p => p.ProteinGroupName).Distinct()) + "\t");
            str.Append(string.Join(";", proteinGroups.Select(p => p.GeneName).Distinct()) + "\t");
            str.Append(string.Join(";", proteinGroups.Select(p => p.Organism).Distinct()) + "\t");

            foreach (var file in rawFiles)
                str.Append(intensities[file] + "\t");
            foreach (var file in rawFiles)
                str.Append(detectionTypes[file] + "\t");

            return str.ToString();
        }

        public override int GetHashCode()
        {
            return Sequence.GetHashCode();
        }

        #endregion Public Methods
    }
}