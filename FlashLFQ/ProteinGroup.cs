using System.Collections.Generic;
using System.Text;

namespace FlashLFQ
{
    public class ProteinGroup
    {
        public readonly string ProteinGroupName;
        public readonly string GeneName;
        public readonly string Organism;
        private Dictionary<SpectraFileInfo, double> Intensities;
        public Dictionary<string, ProteinQuantificationEngineResult> ConditionToQuantificationResults;

        public ProteinGroup(string proteinGroupName, string geneName, string organism)
        {
            ProteinGroupName = proteinGroupName;
            GeneName = geneName;
            Organism = organism;
            Intensities = new Dictionary<SpectraFileInfo, double>();
            ConditionToQuantificationResults = new Dictionary<string, ProteinQuantificationEngineResult>();
        }

        public double GetIntensity(SpectraFileInfo fileInfo)
        {
            if (Intensities.TryGetValue(fileInfo, out double intensity))
            {
                return intensity;
            }
            else
            {
                return 0;
            }
        }

        public void SetIntensity(SpectraFileInfo fileInfo, double intensity)
        {
            if (Intensities.ContainsKey(fileInfo))
            {
                Intensities[fileInfo] = intensity;
            }
            else
            {
                Intensities.Add(fileInfo, intensity);
            }
        }

        public static string TabSeparatedHeader(List<SpectraFileInfo> rawFiles)
        {
            var sb = new StringBuilder();
            sb.Append("Protein Groups" + "\t");
            sb.Append("Gene Name" + "\t");
            sb.Append("Organism" + "\t");

            foreach (var rawfile in rawFiles)
            {
                sb.Append("Intensity_" + rawfile.FilenameWithoutExtension + "\t");
            }

            return sb.ToString();
        }

        public string ToString(List<SpectraFileInfo> rawFiles)
        {
            StringBuilder str = new StringBuilder();
            str.Append(ProteinGroupName + "\t");
            str.Append(GeneName + "\t");
            str.Append(Organism + "\t");

            foreach (var file in rawFiles)
            {
                if (Intensities.TryGetValue(file, out double intensity))
                {
                    str.Append(intensity + "\t");
                }
                else
                {
                    str.Append(0 + "\t");
                }
            }

            return str.ToString();
        }

        public override bool Equals(object obj)
        {
            return ((ProteinGroup)obj).ProteinGroupName.Equals(this.ProteinGroupName);
        }

        public override int GetHashCode()
        {
            return ProteinGroupName.GetHashCode();
        }
    }
}