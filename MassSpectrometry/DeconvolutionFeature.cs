using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    public class DeconvolutionFeature
    {
        #region Public Fields

        public List<(int scanNumber, double elutionTime, IsotopicEnvelope isotopicEnvelope)> isotopicEnvelopes = new List<(int scanNumber, double elutionTime, IsotopicEnvelope isotopicEnvelope)>();

        #endregion Public Fields

        #region Public Properties

        public double Mass { get; private set; }

        public int NumPeaks { get; private set; }

        public IEnumerable<int> AllCharges { get { return isotopicEnvelopes.Select(b => b.isotopicEnvelope.charge).ToList(); } }

        #endregion Public Properties

        #region Internal Methods

        internal void AddEnvelope(int scanNumber, double elutionTime, IsotopicEnvelope isotopicEnvelope)
        {
            isotopicEnvelopes.Add((scanNumber, elutionTime, isotopicEnvelope));
            Mass = isotopicEnvelopes.Select(b => b.isotopicEnvelope.monoisotopicMass).Average();
            NumPeaks = isotopicEnvelopes.Select(b => b.isotopicEnvelope.peaks.Count).Sum();
        }

        #endregion Internal Methods
    }
}