using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FlashLFQ
{
    public class Peptide
    {
        public readonly string Sequence;
        public readonly string BaseSequence;
        public readonly int? PeakOrder; // only used in IsoTracker mode
        public readonly int? IsoGroupIndex; // only used in IsoTracker mode for output grouping
        private Dictionary<SpectraFileInfo, double> Intensities;
        private Dictionary<SpectraFileInfo, double> RetentionTimes;
        private Dictionary<SpectraFileInfo, DetectionType> DetectionTypes;
        public readonly HashSet<ProteinGroup> ProteinGroups;
        public readonly bool UseForProteinQuant;
        public double IonizationEfficiency;

        public Peptide(string sequence, string baseSequence, bool useForProteinQuant, HashSet<ProteinGroup> proteinGroups)
        {
            Sequence = sequence;
            BaseSequence = baseSequence;
            Intensities = new Dictionary<SpectraFileInfo, double>();
            RetentionTimes = new Dictionary<SpectraFileInfo, double>();
            DetectionTypes = new Dictionary<SpectraFileInfo, DetectionType>();
            this.ProteinGroups = proteinGroups;
            this.UseForProteinQuant = useForProteinQuant;
        }

        /// <summary>
        /// The constructor for IsoPeptide
        /// </summary>
        /// <param name="sequence"></param>
        /// <param name="baseSequence"></param>
        /// <param name="useForProteinQuant"></param>
        /// <param name="proteinGroups"></param>
        /// <param name="isoGroupIndex"></param>
        /// <param name="peakOrder"></param>
        public Peptide(string sequence, string baseSequence, bool useForProteinQuant, HashSet<ProteinGroup> proteinGroups, int isoGroupIndex,int peakOrder):
            this(sequence, baseSequence, useForProteinQuant, proteinGroups)
        {
            IsoGroupIndex = isoGroupIndex;
            PeakOrder = peakOrder;
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
            return sb.ToString().TrimEnd('\t');
        }

        public static string TabSeparatedHeader_IsoTracker(List<SpectraFileInfo> rawFiles)
        {
            var sb = new StringBuilder();
            sb.Append("Sequence" + "\t");
            sb.Append("Base Sequence" + "\t");
            sb.Append("Peak Order" + "\t");
            sb.Append("Protein Groups" + "\t");
            sb.Append("Gene Names" + "\t");
            sb.Append("Organism" + "\t");
            foreach (var rawfile in rawFiles)
            {
                sb.Append("Intensity_" + rawfile.FilenameWithoutExtension + "\t");
            }
            foreach (var rawfile in rawFiles)
            {
                sb.Append("RetentionTime (min)_" + rawfile.FilenameWithoutExtension + "\t");
            }
            foreach (var rawfile in rawFiles)
            {
                sb.Append("Detection Type_" + rawfile.FilenameWithoutExtension + "\t");
            }
            return sb.ToString().TrimEnd('\t');
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

        public double GetRetentionTime(SpectraFileInfo fileInfo)
        {
            if (RetentionTimes.TryGetValue(fileInfo, out double retentionTime))
            {
                return retentionTime;
            }
            else
            {
                return 0;
            }
        }

        public void SetRetentionTime(SpectraFileInfo fileInfo, double retentionTime) 
        {
            if (RetentionTimes.ContainsKey(fileInfo))
            {
                RetentionTimes[fileInfo] = retentionTime;
            }
            else
            {
                RetentionTimes.Add(fileInfo, retentionTime);
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

        public void SetIsobaricPeptide(List<ChromatographicPeak> peakList)
        {
            foreach (var peak in peakList.Where(p=>p != null))
            {
                RetentionTimes[peak.SpectraFileInfo] = peak.ApexRetentionTime;
                Intensities[peak.SpectraFileInfo] = peak.Apex.Intensity;
                DetectionTypes[peak.SpectraFileInfo] = peak.DetectionType;
            }
        }



        public string ToString(List<SpectraFileInfo> rawFiles, bool IsoTracker = false)
        {
            if (IsoTracker) // For IsoTracker mode, we add the retention time to the output
            {
                StringBuilder str = new StringBuilder();
                str.Append(Sequence + "\t");
                str.Append(BaseSequence + "\t");
                str.Append(PeakOrder != null ? PeakOrder + "\t" : "" + "\t");

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
                    double intensity = GetIntensity(file);
                    str.Append(intensity + "\t");
                }

                foreach (var file in rawFiles)
                {
                    double Rt = GetRetentionTime(file);
                    str.Append(Rt + "\t");
                }

                foreach (var file in rawFiles)
                {
                    DetectionType detectionType = GetDetectionType(file);
                    str.Append(detectionType + "\t");
                }

                return str.ToString().TrimEnd('\t');
            }

            else 
            {
                StringBuilder str = new StringBuilder();
                str.Append(Sequence + "\t");
                //str.Append("N/A" + "\t"); why we need this?
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

                return str.ToString().TrimEnd('\t');
            }
            
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