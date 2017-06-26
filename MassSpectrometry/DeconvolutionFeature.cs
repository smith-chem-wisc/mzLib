using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MassSpectrometry
{
    internal class DeconvolutionFeature
    {

        #region Private Fields

        private List<IsotopicEnvelope> isotopicEnvelopes = new List<IsotopicEnvelope>();

        #endregion Private Fields

        #region Public Properties

        public double Mass { get { return isotopicEnvelopes.Select(b => b.monoisotopicMass).Average(); } }

        public int NumPeaks { get { return isotopicEnvelopes.Select(b => b.peaks.Count).Sum(); } }

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
        }

        #endregion Internal Methods

    }
}