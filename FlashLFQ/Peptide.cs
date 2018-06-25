using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FlashLFQ
{
    public enum DetectionType { MSMS, MBR, NotDetected, MSMSAmbiguousPeakfinding, MSMSIdentifiedButNotQuantified }

    public class Peptide
    {
        public readonly string Sequence;
        private Dictionary<SpectraFileInfo, double> intensities;
        private Dictionary<SpectraFileInfo, DetectionType> detectionTypes;
        public HashSet<ProteinGroup> proteinGroups;

        public Peptide(string sequence)
        {
            this.Sequence = sequence;
            intensities = new Dictionary<SpectraFileInfo, double>();
            detectionTypes = new Dictionary<SpectraFileInfo, DetectionType>();
            proteinGroups = new HashSet<ProteinGroup>();
        }

        public static string TabSeparatedHeader(List<SpectraFileInfo> rawFiles)
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

        public double GetIntensity(SpectraFileInfo fileInfo)
        {
            if (intensities.TryGetValue(fileInfo, out double intensity))
                return intensity;
            else
                return 0;
        }

        public void SetIntensity(SpectraFileInfo fileInfo, double intensity)
        {
            if (intensities.ContainsKey(fileInfo))
                intensities[fileInfo] = intensity;
            else
                intensities.Add(fileInfo, intensity);
        }

        public DetectionType GetDetectionType(SpectraFileInfo fileInfo)
        {
            if (detectionTypes.TryGetValue(fileInfo, out DetectionType detectionType))
                return detectionType;
            else
                return DetectionType.NotDetected;
        }

        public void SetDetectionType(SpectraFileInfo fileInfo, DetectionType detectionType)
        {
            if (detectionTypes.ContainsKey(fileInfo))
                detectionTypes[fileInfo] = detectionType;
            else
                detectionTypes.Add(fileInfo, detectionType);
        }

        public string ToString(List<SpectraFileInfo> rawFiles)
        {
            StringBuilder str = new StringBuilder();
            str.Append(Sequence + "\t");
            str.Append(string.Join(";", proteinGroups.Select(p => p.ProteinGroupName).Distinct()) + "\t");
            str.Append(string.Join(";", proteinGroups.Select(p => p.GeneName).Distinct()) + "\t");
            str.Append(string.Join(";", proteinGroups.Select(p => p.Organism).Distinct()) + "\t");

            foreach (var file in rawFiles)
            {
                str.Append(GetIntensity(file) + "\t");
            }
            foreach (var file in rawFiles)
            {
                str.Append(GetDetectionType(file) + "\t");
            }

            return str.ToString();
        }

        public override bool Equals(object obj)
        {
            return ((Peptide)obj).Sequence.Equals(this.Sequence);
        }

        public override int GetHashCode()
        {
            return Sequence.GetHashCode();
        }
    }
}