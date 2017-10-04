﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MassSpectrometry
{
    public class DeconvolutionFeatureWithMassesAndScans
    {
        #region Private Fields

        private List<DeconvolutionFeature> groups = new List<DeconvolutionFeature>();

        #endregion Private Fields

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

        public int NumPeaks
        {
            get { return groups.Select(b => b.NumPeaks).Sum(); }
        }

        public double MinElutionTime { get; private set; }
        public double MaxElutionTime { get; private set; }
        public double TotalIntensity { get; private set; }

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
                + NumPeaks + "\t"
                + (MaxScanIndex - MinScanIndex + 1) + "\t"
                + MinScanIndex + "\t"
                + MaxScanIndex + "\t"
                + ((MinElutionTime + MaxElutionTime) / 2).ToString("F2") + "\t"
                + TotalIntensity.ToString("E5") + "\t"
                + string.Join(",", new HashSet<int>(groups.SelectMany(b => b.AllCharges)).OrderBy(b => b)) + "\t"
                + groups.OrderByDescending(b => b.MostIntenseEnvelope.totalIntensity).First().MostIntenseEnvelope.ToString();
        }

        public void AddEnvelope(IsotopicEnvelope isotopicEnvelope, int scanIndex, double elutionTime)
        {
            MinScanIndex = Math.Min(scanIndex, MinScanIndex);
            MaxScanIndex = Math.Max(scanIndex, MaxScanIndex);
            MinElutionTime = Math.Min(elutionTime, MinElutionTime);
            MaxElutionTime = Math.Max(elutionTime, MaxElutionTime);
            foreach (var massGroup in groups)
            {
                if (Math.Abs(massGroup.Mass - isotopicEnvelope.monoisotopicMass) < 0.5)
                {
                    massGroup.AddEnvelope(isotopicEnvelope);
                    Mass = groups.OrderBy(b => -b.NumPeaks).First().Mass;
                    return;
                }
            }
            var newMassGroup = new DeconvolutionFeature();
            newMassGroup.AddEnvelope(isotopicEnvelope);
            groups.Add(newMassGroup);

            Mass = groups.OrderBy(b => -b.NumPeaks).First().Mass;
            TotalIntensity += isotopicEnvelope.peaks.Sum(b => b.Item2);
        }

        #endregion Public Methods
    }
}