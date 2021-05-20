using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MassSpectrometry
{
    public class DeconvolutionFeatureWithMassesAndScans
    {
        public List<DeconvolutionFeature> groups = new List<DeconvolutionFeature>();

        public DeconvolutionFeatureWithMassesAndScans()
        {
            MinScanIndex = int.MaxValue;

            MaxScanIndex = int.MinValue;
            MinElutionTime = double.MaxValue;

            MaxElutionTime = double.MinValue;
        }

        public int MinScanIndex { get; private set; }
        public int MaxScanIndex { get; private set; }
        public double Mass { get; private set; }

        public double Score
        {
            get
            {
                return Math.Log(
                          Math.Pow(TotalNormalizedIntensity, 0.1)
                        * Math.Pow(Math.Max((MaxElutionTime - MinElutionTime * 60), 1), 0.1)
                        * Math.Pow((new HashSet<int>(groups.SelectMany(b => b.AllCharges)).OrderBy(b => b)).Count(), 1)
                        * Math.Pow((double)groups.Select(b => b.NumPeaks).Sum() / (MaxScanIndex - MinScanIndex + 1), 1));
            }
        }

        public int NumPeaks
        {
            get { return groups.Select(b => b.NumPeaks).Sum(); }
        }

        public double MinElutionTime { get; private set; }
        public double MaxElutionTime { get; private set; }
        public double TotalNormalizedIntensity { get; private set; }
        public IsotopicEnvelope MostIntenseEnvelope { get; private set; }
        public double MostIntenseEnvelopeElutionTime { get; private set; }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append(OneLineString());
            foreach (var heh in groups.OrderBy(b => -b.NumPeaks))
            {
                sb.AppendLine();
                sb.Append("  " + heh.ToString());
            }

            return sb.ToString();
        }

        public string OneLineString()
        {
            List<(double elutionTime, double intensity)> elutionTimeAndIntensity = new List<(double elutionTime, double intensity)>();
            foreach (var (scanNumber, elutionTime, isotopicEnvelope) in groups.SelectMany(b => b.isotopicEnvelopes).Where(b => b.isotopicEnvelope.Charge == MostIntenseEnvelope.Charge))
                elutionTimeAndIntensity.Add((elutionTime, isotopicEnvelope.TotalIntensity));

            int maxCharge = groups.SelectMany(p => p.AllCharges).Max();
            var t = groups.SelectMany(p => p.isotopicEnvelopes);
            string elutionString = "";
            for (int z = 1; z <= maxCharge; z++)
            {
                string str = "[" + z + "|";
                var isotopicEnvelopes = t.Where(p => p.isotopicEnvelope.Charge == z);
                foreach (var (scanNumber, elutionTime, isotopicEnvelope) in isotopicEnvelopes)
                {
                    str += Math.Round(elutionTime, 2) + ";" + isotopicEnvelope.TotalIntensity + ",";
                }
                str += "]";

                elutionString += str;
            }

            var elutionOfMostIntenseCharge = string.Join(";", elutionTimeAndIntensity.OrderBy(b => b.elutionTime).Select(b => b.intensity));

            return Mass.ToString("G8") + "\t"
                + Score + "\t"
                + NumPeaks + "\t"
                + (MaxScanIndex - MinScanIndex + 1) + "\t"
                + MinScanIndex + "\t"
                + MaxScanIndex + "\t"
                + MinElutionTime + "\t"
                + MaxElutionTime + "\t"
                + TotalNormalizedIntensity.ToString("E5") + "\t"
                + string.Join(",", new HashSet<int>(groups.SelectMany(b => b.AllCharges)).OrderBy(b => b)) + "\t"
                + (MostIntenseEnvelopeElutionTime).ToString("F2") + "\t"
                + MostIntenseEnvelope.ToString() + '\t'
                + elutionOfMostIntenseCharge + '\t'
                + elutionString;
        }

        public void AddEnvelope(IsotopicEnvelope isotopicEnvelope, int scanIndex, double elutionTime)
        {
            MinScanIndex = Math.Min(scanIndex, MinScanIndex);
            MaxScanIndex = Math.Max(scanIndex, MaxScanIndex);
            MinElutionTime = Math.Min(elutionTime, MinElutionTime);
            MaxElutionTime = Math.Max(elutionTime, MaxElutionTime);
            bool added = false;
            foreach (var massGroup in groups)
            {
                if (Math.Abs(massGroup.Mass - isotopicEnvelope.MonoisotopicMass) < 0.5)
                {
                    massGroup.AddEnvelope(scanIndex, elutionTime, isotopicEnvelope);
                    added = true;
                    break;
                }
            }
            if (!added)
            {
                var newMassGroup = new DeconvolutionFeature();
                newMassGroup.AddEnvelope(scanIndex, elutionTime, isotopicEnvelope);
                groups.Add(newMassGroup);
            }

            Mass = groups.OrderBy(b => -b.NumPeaks).First().Mass;
            TotalNormalizedIntensity += isotopicEnvelope.TotalIntensity / isotopicEnvelope.Charge;

            if (MostIntenseEnvelope == null || MostIntenseEnvelope.TotalIntensity < isotopicEnvelope.TotalIntensity)
            {
                MostIntenseEnvelope = isotopicEnvelope;
                MostIntenseEnvelopeElutionTime = elutionTime;
            }
        }
    }
}