﻿using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MassSpectrometry
{
    public class DeconvolutionFeature
    {
        #region Public Fields

        public List<IsotopicEnvelope> isotopicEnvelopes = new List<IsotopicEnvelope>();

        #endregion Public Fields

        #region Public Properties

        public double Mass { get; private set; }

        public int NumPeaks { get; private set; }

        public IsotopicEnvelope MostIntenseEnvelope { get { return isotopicEnvelopes.OrderByDescending(b => b.totalIntensity).First(); } }

        public IEnumerable<int> AllCharges { get { return isotopicEnvelopes.Select(b => b.charge).ToList(); } }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append(" Mass: " + Mass + " NumPeaks: " + NumPeaks);
            foreach (var heh in isotopicEnvelopes.OrderBy(b => -b.peaks.Count))
            {
                sb.AppendLine();
                sb.Append("     " + heh.ToString());
            }
            return sb.ToString();
        }

        #endregion Public Methods

        #region Internal Methods

        internal void AddEnvelope(IsotopicEnvelope isotopicEnvelope)
        {
            isotopicEnvelopes.Add(isotopicEnvelope);
            Mass = isotopicEnvelopes.Select(b => b.monoisotopicMass).Average();
            NumPeaks = isotopicEnvelopes.Select(b => b.peaks.Count).Sum();
        }

        #endregion Internal Methods
    }
}