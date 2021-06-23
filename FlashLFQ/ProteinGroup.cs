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
        public List<Peptide> PeptidesUsedForProteinQuantification { get; private set; }

        public ProteinGroup(string proteinGroupName, string geneName, string organism)
        {
            ProteinGroupName = proteinGroupName;
            GeneName = geneName;
            Organism = organism;
            Intensities = new Dictionary<SpectraFileInfo, double>();
            ConditionToQuantificationResults = new Dictionary<string, ProteinQuantificationEngineResult>();
            PeptidesUsedForProteinQuantification = new List<Peptide>();
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

        public (double intensity, DetectionType detType)[,] GetPeptideIntensityAndDetectionTypeGrid(List<Peptide> peptides, out List<(string, int)> bioreps)
        {
            PeptidesUsedForProteinQuantification = peptides;
            bioreps = new List<(string, int)>();

            foreach (var condition in Intensities.Keys.GroupBy(p => p.Condition).OrderBy(p => p.Key))
            {
                foreach (int biorep in condition.OrderBy(p => p.BiologicalReplicate).Select(p => p.BiologicalReplicate).Distinct())
                {
                    bioreps.Add((condition.Key, biorep));
                }
            }

            var grid = new (double, DetectionType)[peptides.Count, bioreps.Count];

            int b = 0;
            foreach (var condition in Intensities.Keys.GroupBy(p => p.Condition).OrderBy(p => p.Key))
            {
                foreach (var biorep in condition.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                {
                    for (int pep = 0; pep < peptides.Count; pep++)
                    {
                        Peptide peptide = peptides[pep];

                        List<(double intensity, DetectionType detType)> fractionSummaries = new List<(double, DetectionType)>();

                        foreach (var fraction in biorep.GroupBy(p => p.Fraction))
                        {
                            DetectionType bestDetectionType = DetectionType.NotDetected;
                            var techrepMeasurements = new List<(double intensity, DetectionType detType)> { (0, DetectionType.NotDetected) };

                            foreach (SpectraFileInfo replicate in fraction.OrderBy(p => p.TechnicalReplicate))
                            {
                                double replicateIntensity = peptide.GetIntensity(replicate);
                                DetectionType replicateDetectionType = peptide.GetDetectionType(replicate);

                                if (replicateIntensity > 0)
                                {
                                    techrepMeasurements.RemoveAll(p => p.detType == DetectionType.NotDetected);

                                    if (replicateDetectionType == DetectionType.MSMS)
                                    {
                                        techrepMeasurements.RemoveAll(p => p.detType != DetectionType.MSMS);
                                        techrepMeasurements.Add((replicateIntensity, replicateDetectionType));
                                        bestDetectionType = DetectionType.MSMS;
                                    }
                                    else if (replicateDetectionType == DetectionType.MBR && bestDetectionType != DetectionType.MSMS)
                                    {
                                        techrepMeasurements.RemoveAll(p => p.detType == DetectionType.Imputed);
                                        techrepMeasurements.Add((replicateIntensity, replicateDetectionType));
                                        bestDetectionType = DetectionType.MBR;
                                    }
                                    else if (replicateDetectionType == DetectionType.Imputed &&
                                        (bestDetectionType == DetectionType.NotDetected || bestDetectionType == DetectionType.Imputed))
                                    {
                                        techrepMeasurements.Add((replicateIntensity, replicateDetectionType));
                                        bestDetectionType = DetectionType.Imputed;
                                    }
                                }
                            }

                            (double intensity, DetectionType detType) fractionSummary =
                                (techrepMeasurements.Where(p => p.detType == bestDetectionType).Average(p => p.intensity), bestDetectionType);

                            fractionSummaries.Add(fractionSummary);
                        }

                        (double intensity, DetectionType detType) bestFractionSummary = fractionSummaries.OrderBy(p => p.detType).ThenByDescending(p => p.intensity).First();
                        grid[pep, b] = bestFractionSummary;
                    }

                    b++;
                }
            }

            return grid;
        }

        public static string TabSeparatedHeader(List<SpectraFileInfo> spectraFiles)
        {
            var sb = new StringBuilder();
            sb.Append("Protein Groups" + "\t");
            sb.Append("Gene Name" + "\t");
            sb.Append("Organism" + "\t");
            sb.Append("Number Quantified Peptides" + "\t");

            bool unfractionated = spectraFiles.Select(p => p.Fraction).Distinct().Count() == 1;
            bool conditionsUndefined = spectraFiles.All(p => p.Condition == "Default") || spectraFiles.All(p => string.IsNullOrWhiteSpace(p.Condition));

            foreach (var sampleGroup in spectraFiles.GroupBy(p => p.Condition))
            {
                foreach (var sample in sampleGroup.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                {
                    if (conditionsUndefined && unfractionated)
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
            sb.Append(PeptidesUsedForProteinQuantification.Count + "\t");

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

        public string ToPeptideIntensityAndDetectionTypeGridString()
        {
            StringBuilder sb = new StringBuilder();
            var grid = GetPeptideIntensityAndDetectionTypeGrid(PeptidesUsedForProteinQuantification, out var bioreps);

            var header = 
                "\t" + 
                string.Join('\t', bioreps.Select(p => p.Item1 + "_" + (p.Item2 + 1) + "_Intensity")) +
                "\t" +
                string.Join('\t', bioreps.Select(p => p.Item1 + "_" + (p.Item2 + 1) + "_DetectionType"));
            sb.Append(header);
            sb.Append('\n');

            for (int i = 0; i < grid.GetLength(0); i++)
            {
                sb.Append(PeptidesUsedForProteinQuantification[i].Sequence);
                sb.Append('\t');

                for (int j = 0; j < grid.GetLength(1); j++)
                {
                    sb.Append(grid[i, j].intensity);
                    sb.Append('\t');
                }

                for (int j = 0; j < grid.GetLength(1); j++)
                {
                    sb.Append(grid[i, j].detType);
                    sb.Append('\t');
                }

                sb.Append('\n');
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