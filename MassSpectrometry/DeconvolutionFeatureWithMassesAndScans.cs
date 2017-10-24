using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MassSpectrometry
{
    public class DeconvolutionFeatureWithMassesAndScans
    {
        #region Public Fields

        public List<DeconvolutionFeature> groups = new List<DeconvolutionFeature>();

        #endregion Public Fields

        #region Public Constructors

        public DeconvolutionFeatureWithMassesAndScans()
        {
            MinScanIndex = int.MaxValue;

            MaxScanIndex = int.MinValue;
            MinElutionTime = double.MaxValue;

            MaxElutionTime = double.MinValue;
        }

        #endregion Public Constructors

        #region Public Properties

        public int MinScanIndex { get; private set; }

        public int MaxScanIndex { get; private set; }

        public double Mass { get; private set; }

        public double Score
        {
            get
            {
                return Math.Log(
                          Math.Pow(TotalIntensity, 0.1)
                        * Math.Pow(Math.Max((MaxElutionTime - MinElutionTime * 60), 1), 0.1)
                        * Math.Pow((new HashSet<int>(groups.SelectMany(b => b.AllCharges)).OrderBy(b => b)).Count(), 2));
            }
        }

        public int NumPeaks
        {
            get { return groups.Select(b => b.NumPeaks).Sum(); }
        }

        public double MinElutionTime { get; private set; }
        public double MaxElutionTime { get; private set; }
        public double TotalIntensity { get; private set; }
        public IsotopicEnvelope MostIntenseEnvelope { get; private set; }
        public double MostIntenseEnvelopeElutionTime { get; private set; }

        #endregion Public Properties

        #region Public Methods

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
            return Mass.ToString("G8") + "\t"
                + Score + "\t"
                + NumPeaks + "\t"
                + (MaxScanIndex - MinScanIndex + 1) + "\t"
                + MinScanIndex + "\t"
                + MaxScanIndex + "\t"
                + TotalIntensity.ToString("E5") + "\t"
                + string.Join(",", new HashSet<int>(groups.SelectMany(b => b.AllCharges)).OrderBy(b => b)) + "\t"
                + (MostIntenseEnvelopeElutionTime).ToString("F2") + "\t"
                + groups.OrderByDescending(b => b.MostIntenseEnvelope.totalIntensity).First().MostIntenseEnvelope.ToString();
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
                if (Math.Abs(massGroup.Mass - isotopicEnvelope.monoisotopicMass) < 0.5)
                {
                    massGroup.AddEnvelope(isotopicEnvelope);
                    added = true;
                    break;
                }
            }
            if (!added)
            {
                var newMassGroup = new DeconvolutionFeature();
                newMassGroup.AddEnvelope(isotopicEnvelope);
                groups.Add(newMassGroup);
            }

            Mass = groups.OrderBy(b => -b.NumPeaks).First().Mass;
            TotalIntensity += isotopicEnvelope.peaks.Sum(b => b.Item2);

            if (MostIntenseEnvelope.totalIntensity < isotopicEnvelope.totalIntensity)
            {
                MostIntenseEnvelope = isotopicEnvelope;
                MostIntenseEnvelopeElutionTime = elutionTime;
            }
        }

        #endregion Public Methods
    }
}