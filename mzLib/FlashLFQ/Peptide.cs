using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FlashLFQ
{
    public class Peptide
    {
        public readonly string Sequence;
        public readonly string BaseSequence;
        private Dictionary<SpectraFileInfo, double> Intensities;
        private Dictionary<SpectraFileInfo, DetectionType> DetectionTypes;
        public readonly HashSet<ProteinGroup> ProteinGroups;
        public readonly bool UseForProteinQuant;
        public double IonizationEfficiency;

        public Peptide(string sequence, string baseSequence, bool useForProteinQuant, HashSet<ProteinGroup> proteinGroups)
        {
            Sequence = sequence;
            BaseSequence = baseSequence;
            Intensities = new Dictionary<SpectraFileInfo, double>();
            DetectionTypes = new Dictionary<SpectraFileInfo, DetectionType>();
            this.ProteinGroups = proteinGroups;
            this.UseForProteinQuant = useForProteinQuant;
        }

        public static string TabSeparatedHeader(List<SpectraFileInfo> rawFiles)
        {
            var sb = new StringBuilder();
            sb.Append("Sequence" + "\t");
            sb.Append("Base Sequence" + "\t");
            sb.Append("Protein Groups" + "\t");
            sb.Append("Gene Names" + "\t");
            sb.Append("Organism" + "\t");
            foreach (var rawfile in rawFiles)
            {
                sb.Append("Intensity_" + rawfile.FilenameWithoutExtension + "\t");
            }
            foreach (var rawfile in rawFiles)
            {
                sb.Append("Detection Type_" + rawfile.FilenameWithoutExtension + "\t");
            }
            return sb.ToString();
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

        public DetectionType GetDetectionType(SpectraFileInfo fileInfo)
        {
            if (DetectionTypes.TryGetValue(fileInfo, out DetectionType detectionType))
            {
                return detectionType;
            }
            else
            {
                return DetectionType.NotDetected;
            }
        }

        public void SetDetectionType(SpectraFileInfo fileInfo, DetectionType detectionType)
        {
            if (DetectionTypes.ContainsKey(fileInfo))
            {
                DetectionTypes[fileInfo] = detectionType;
            }
            else
            {
                DetectionTypes.Add(fileInfo, detectionType);
            }
        }

        public string ToString(List<SpectraFileInfo> rawFiles)
        {
            StringBuilder str = new StringBuilder();
            str.Append(Sequence + "\t");
            str.Append(BaseSequence + "\t");

            var orderedProteinGroups = ProteinGroups.OrderBy(p => p.ProteinGroupName).ToList();

            var proteinsCount = orderedProteinGroups.Select(p => p.ProteinGroupName).Distinct().Count();
            var genesCount = orderedProteinGroups.Select(p => p.GeneName).Distinct().Count();
            var organismsCount = orderedProteinGroups.Select(p => p.Organism).Distinct().Count();

            str.Append(proteinsCount > 1 ? string.Join(";", orderedProteinGroups.Select(p => p.ProteinGroupName)) + "\t" :
                orderedProteinGroups.Any() ? orderedProteinGroups.First().ProteinGroupName + "\t" : "\t");

            str.Append(genesCount > 1 ? string.Join(";", orderedProteinGroups.Select(p => p.GeneName)) + "\t" :
                orderedProteinGroups.Any() ? orderedProteinGroups.First().GeneName + "\t" : "\t");

            str.Append(organismsCount > 1 ? string.Join(";", orderedProteinGroups.Select(p => p.Organism)) + "\t" :
                orderedProteinGroups.Any() ? orderedProteinGroups.First().Organism + "\t" : "\t");

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

        public bool UnambiguousPeptideQuant()
        {
            return Intensities.Values.Any(p => p > 0) && DetectionTypes.Values.Any(p => p != DetectionType.MSMSAmbiguousPeakfinding);
        }
    }
}