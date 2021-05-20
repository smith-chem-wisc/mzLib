using System.Collections.Generic;
using System.Linq;
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

        public static string TabSeparatedHeader(List<SpectraFileInfo> spectraFiles)
        {
            var sb = new StringBuilder();
            sb.Append("Protein Groups" + "\t");
            sb.Append("Gene Name" + "\t");
            sb.Append("Organism" + "\t");

            bool unfractionated = spectraFiles.Select(p => p.Fraction).Distinct().Count() == 1;
            bool conditionsDefined = spectraFiles.All(p => p.Condition == "Default") || spectraFiles.All(p => string.IsNullOrWhiteSpace(p.Condition));

            foreach (var sampleGroup in spectraFiles.GroupBy(p => p.Condition))
            {
                foreach (var sample in sampleGroup.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                {
                    if (!conditionsDefined && unfractionated)
                    {
                        sb.Append("Intensity_" + sample.First().FilenameWithoutExtension + "\t");
                    }
                    else
                    {
                        sb.Append("Intensity_" + sample.First().Condition + "_" + (sample.First().BiologicalReplicate + 1) + "\t");
                    }
                }
            }

            return sb.ToString();
        }

        public string ToString(List<SpectraFileInfo> spectraFiles)
        {
            StringBuilder sb = new StringBuilder();
            sb.Append(ProteinGroupName + "\t");
            sb.Append(GeneName + "\t");
            sb.Append(Organism + "\t");

            bool unfractionated = spectraFiles.Select(p => p.Fraction).Distinct().Count() == 1;
            bool conditionsDefined = spectraFiles.All(p => p.Condition == "Default") || spectraFiles.All(p => string.IsNullOrWhiteSpace(p.Condition));

            foreach (var sampleGroup in spectraFiles.GroupBy(p => p.Condition))
            {
                foreach (var sample in sampleGroup.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                {
                    if (!conditionsDefined && unfractionated)
                    {
                        sb.Append(GetIntensity(sample.First()) + "\t");
                    }
                    else
                    {
                        double summedIntensity = sample.Sum(p => GetIntensity(p));
                        sb.Append(summedIntensity + "\t");
                    }
                }
            }

            return sb.ToString();
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