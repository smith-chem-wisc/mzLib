using System.Collections.Generic;
using System.Linq;
using System.Text;
using MassSpectrometry;

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

        /// <summary>
        /// Whether a sample's column should be labelled with its file name rather than its
        /// condition and biological replicate.
        /// </summary>
        /// <remarks>
        /// A file name is the only meaningful label when no experimental design was supplied — every
        /// condition is blank or "Default", so a condition label would read "Intensity__1". Once real
        /// conditions exist, the condition and replicate identify the sample and the file name does
        /// not. A file name is also meaningless for fractionated data, where one sample spans several
        /// files and <c>sample.First()</c> names only one of them.
        /// <para>
        /// The header and the row must agree on this, so both call here rather than each recomputing
        /// it.
        /// </para>
        /// </remarks>
        private static bool LabelSamplesByFileName(List<SpectraFileInfo> spectraFiles)
        {
            bool unfractionated = spectraFiles.Select(p => p.Fraction).Distinct().Count() == 1;
            bool conditionsUndefined = spectraFiles.All(p => p.Condition == "Default") || spectraFiles.All(p => string.IsNullOrWhiteSpace(p.Condition));

            return conditionsUndefined && unfractionated;
        }

        public static string TabSeparatedHeader(List<SpectraFileInfo> spectraFiles)
        {
            var sb = new StringBuilder();
            sb.Append("Protein Groups" + "\t");
            sb.Append("Gene Name" + "\t");
            sb.Append("Organism" + "\t");

            bool labelByFileName = LabelSamplesByFileName(spectraFiles);

            foreach (var sampleGroup in spectraFiles.GroupBy(p => p.Condition))
            {
                foreach (var sample in sampleGroup.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                {
                    if (labelByFileName)
                    {
                        sb.Append("Intensity_" + sample.First().FilenameWithoutExtension + "\t");
                    }
                    else
                    {
                        sb.Append("Intensity_" + sample.First().Condition + "_" + (sample.First().BiologicalReplicate + 1) + "\t");
                    }
                }
            }

            return sb.ToString().TrimEnd('\t');
        }

        public string ToString(List<SpectraFileInfo> spectraFiles)
        {
            StringBuilder sb = new StringBuilder();
            sb.Append(ProteinGroupName + "\t");
            sb.Append(GeneName + "\t");
            sb.Append(Organism + "\t");

            bool labelByFileName = LabelSamplesByFileName(spectraFiles);

            foreach (var sampleGroup in spectraFiles.GroupBy(p => p.Condition))
            {
                foreach (var sample in sampleGroup.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                {
                    if (labelByFileName)
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
            
            return sb.ToString().TrimEnd('\t');
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