using System.Collections.Generic;

namespace MassSpectrometry
{
    public class IsotopicEnvelope
    {
        #region Public Fields

        public readonly List<(double, double)> peaks;
        public readonly double monoisotopicMass;
        public readonly int charge;
        public readonly double totalIntensity;
        public readonly double stDev;
        public readonly int massIndex;

        #endregion Public Fields

        #region Public Constructors

        public IsotopicEnvelope(List<(double, double)> bestListOfPeaks, double bestMonoisotopicMass, int bestChargeState, double bestTotalIntensity, double bestStDev, int bestMassIndex)
        {
            this.peaks = bestListOfPeaks;
            this.monoisotopicMass = bestMonoisotopicMass;
            this.charge = bestChargeState;
            this.totalIntensity = bestTotalIntensity;
            this.stDev = bestStDev;
            this.massIndex = bestMassIndex;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            return charge + "\t" + peaks[0].Item1.ToString("G8") + "\t" + peaks.Count + "\t" + totalIntensity;
        }

        #endregion Public Methods
    }
}