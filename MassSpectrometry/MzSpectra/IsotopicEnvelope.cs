using System.Collections.Generic;

namespace MassSpectrometry
{
    public class IsotopicEnvelope
    {
        public readonly List<(double mz, double intensity)> peaks;
        public readonly double monoisotopicMass;
        public readonly int charge;
        public readonly double totalIntensity;
        public readonly double stDev;
        public readonly int massIndex;

        public IsotopicEnvelope(List<(double mz, double intensity)> bestListOfPeaks, double bestMonoisotopicMass, int bestChargeState, double bestTotalIntensity, double bestStDev, int bestMassIndex)
        {
            this.peaks = bestListOfPeaks;
            this.monoisotopicMass = bestMonoisotopicMass;
            this.charge = bestChargeState;
            this.totalIntensity = bestTotalIntensity;
            this.stDev = bestStDev;
            this.massIndex = bestMassIndex;
        }

        public override string ToString()
        {
            return charge + "\t" + peaks[0].mz.ToString("G8") + "\t" + peaks.Count + "\t" + totalIntensity;
        }
    }
}