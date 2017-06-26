using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Deconvolution
{
    internal class AggregatedIsotopeEnvelopes
    {

        #region Public Fields

        public int minScanIndex = int.MaxValue;

        public int maxScanIndex = int.MinValue;

        public double mass;

        #endregion Public Fields

        #region Private Fields

        private List<AggregatedIsotopeEnvelopesWithSameMass> groups = new List<AggregatedIsotopeEnvelopesWithSameMass>();

        #endregion Private Fields

        #region Public Properties

        public int NumPeaks
        {
            get { return groups.Select(b => b.NumPeaks).Sum(); }
        }

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

        #endregion Public Methods

        #region Internal Methods

        internal string OneLineString()
        {
            return "Mass: " + mass + " NumPeaks: " + NumPeaks + " NumScans: " + (maxScanIndex - minScanIndex + 1) + " ScanRange: " + minScanIndex + " to " + maxScanIndex;
        }

        internal void AddEnvelope(IsotopicEnvelope isotopicEnvelope, int scanIndex)
        {
            minScanIndex = Math.Min(scanIndex, minScanIndex);
            maxScanIndex = Math.Max(scanIndex, maxScanIndex);
            foreach (var massGroup in groups)
            {
                if (Math.Abs(massGroup.Mass - isotopicEnvelope.monoisotopicMass) < 0.5)
                {
                    massGroup.AddEnvelope(isotopicEnvelope);
                    mass = groups.OrderBy(b => -b.NumPeaks).First().Mass;
                    return;
                }
            }
            var newMassGroup = new AggregatedIsotopeEnvelopesWithSameMass();
            newMassGroup.AddEnvelope(isotopicEnvelope);
            groups.Add(newMassGroup);

            mass = groups.OrderBy(b => -b.NumPeaks).First().Mass;
        }

        #endregion Internal Methods

    }
}