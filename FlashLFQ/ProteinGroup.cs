using System.Collections.Generic;
using System.Text;

namespace FlashLFQ
{
    public class ProteinGroup
    {
        public readonly string ProteinGroupName;
        public readonly string GeneName;
        public readonly string Organism;
        public readonly Dictionary<SpectraFileInfo, double> intensities;

        public ProteinGroup(string proteinGroupName, string GeneName, string Organism)
        {
            this.ProteinGroupName = proteinGroupName;
            this.GeneName = GeneName;
            this.Organism = Organism;
            intensities = new Dictionary<SpectraFileInfo, double>();
        }

        public static string TabSeparatedHeader(List<SpectraFileInfo> rawFiles)
        {
            var sb = new StringBuilder();
            sb.Append("Protein Groups" + "\t");
            sb.Append("Gene Name" + "\t");
            sb.Append("Organism" + "\t");

            foreach (var rawfile in rawFiles)
                sb.Append("Intensity_" + rawfile.filenameWithoutExtension + "\t");

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
                if (intensities.TryGetValue(file, out double intensity))
                    str.Append(intensity + "\t");
                else
                    str.Append(0 + "\t");
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