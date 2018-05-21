using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    internal class GeneratedPeak : MzPeakZR
    {
        #region Private Fields

        private List<double> mzs = new List<double>();
        private List<double> intensities = new List<double>();

        #endregion Private Fields

        #region Public Constructors

        public GeneratedPeak(double Mz, double Intensity) : base(Mz, Intensity)
        {
            mzs.Add(Mz);
            intensities.Add(Intensity);
        }

        #endregion Public Constructors

        #region Internal Methods

        internal void AddMzPeak(double anotherMz, double anotherIntensity)
        {
            mzs.Add(anotherMz);
            intensities.Add(anotherIntensity);
            Intensity = intensities.Sum();
            double weightedSumMz = 0;
            for (int i = 0; i < mzs.Count; i++)
                weightedSumMz += mzs[i] * intensities[i];
            Mz = weightedSumMz / Intensity;
        }

        #endregion Internal Methods
    }
}