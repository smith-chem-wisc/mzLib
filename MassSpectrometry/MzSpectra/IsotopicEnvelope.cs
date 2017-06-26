using System;
using System.Collections.Generic;

namespace MassSpectrometry
{
    public class IsotopicEnvelope
    {

        #region Public Fields

        public readonly List<IMzPeak> peaks;
        public readonly int charge;
        public readonly int indexOfMostIntenseHere;
        public readonly double monoisotopicMass;
        public readonly List<double> ratios;

        #endregion Public Fields

        #region Public Constructors

        public IsotopicEnvelope(Tuple<List<IMzPeak>, int, int, double, List<double>> longest)
        {
            this.peaks = longest.Item1;
            this.charge = longest.Item2;
            this.indexOfMostIntenseHere = longest.Item3;
            this.monoisotopicMass = longest.Item4;
            this.ratios = longest.Item5;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            return "MM: " + monoisotopicMass + " charge: " + charge + " numPeaks: " + peaks.Count + " mostIntensePeak: " + peaks[0].Mz;
        }

        #endregion Public Methods

    }
}