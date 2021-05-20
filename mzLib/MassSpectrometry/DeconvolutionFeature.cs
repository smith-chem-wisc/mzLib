using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    public class DeconvolutionFeature
    {
        public List<(int scanNumber, double elutionTime, IsotopicEnvelope isotopicEnvelope)> isotopicEnvelopes = new List<(int scanNumber, double elutionTime, IsotopicEnvelope isotopicEnvelope)>();

        public double Mass { get; private set; }

        public int NumPeaks { get; private set; }

        public IEnumerable<int> AllCharges { get { return isotopicEnvelopes.Select(b => b.isotopicEnvelope.Charge).ToList(); } }

        internal void AddEnvelope(int scanNumber, double elutionTime, IsotopicEnvelope isotopicEnvelope)
        {
            isotopicEnvelopes.Add((scanNumber, elutionTime, isotopicEnvelope));
            Mass = isotopicEnvelopes.Select(b => b.isotopicEnvelope.MonoisotopicMass).Average();
            NumPeaks = isotopicEnvelopes.Select(b => b.isotopicEnvelope.Peaks.Count).Sum();
        }
    }
}