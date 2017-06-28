using System.Collections.Generic;

namespace MassSpectrometry
{
    public class IsotopicEnvelope
    {

        #region Public Fields

        public readonly List<IMzPeak> listOfPeaks;
        public readonly double monoisotopicMass;
        public readonly int chargeState;
        public readonly double totalIntensity;
        public readonly double stDev;
        public readonly int massIndex;
        public readonly int bestJ;

        #endregion Public Fields

        #region Public Constructors

        public IsotopicEnvelope(List<IMzPeak> bestListOfPeaks, double bestMonoisotopicMass, int bestChargeState, double bestTotalIntensity, double bestStDev, int bestMassIndex, int bestJ)
        {
            this.listOfPeaks = bestListOfPeaks;
            this.monoisotopicMass = bestMonoisotopicMass;
            this.chargeState = bestChargeState;
            this.totalIntensity = bestTotalIntensity;
            this.stDev = bestStDev;
            this.massIndex = bestMassIndex;
            this.bestJ = bestJ;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            return "MM: " + monoisotopicMass + " charge: " + chargeState + " numPeaks: " + listOfPeaks.Count + " mostIntensePeak: " + listOfPeaks[0].Mz + " bestJ: " + bestJ;
        }

        #endregion Public Methods

    }
}