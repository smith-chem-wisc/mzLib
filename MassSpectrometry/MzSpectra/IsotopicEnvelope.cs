﻿using System.Collections.Generic;

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
        public readonly int bestJ;

        #endregion Public Fields

        #region Public Constructors

        public IsotopicEnvelope(List<(double, double)> bestListOfPeaks, double bestMonoisotopicMass, int bestChargeState, double bestTotalIntensity, double bestStDev, int bestMassIndex, int bestJ)
        {
            this.peaks = bestListOfPeaks;
            this.monoisotopicMass = bestMonoisotopicMass;
            this.charge = bestChargeState;
            this.totalIntensity = bestTotalIntensity;
            this.stDev = bestStDev;
            this.massIndex = bestMassIndex;
            this.bestJ = bestJ;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            return charge + "\t" + peaks[0].Item1.ToString("G8") + "\t" + peaks.Count;
        }

        #endregion Public Methods
    }
}