using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    internal class GeneratedPeak : MzPeak
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

        internal void Add(double nextMz, double v)
        {
            mzs.Add(nextMz);
            intensities.Add(v);
            Y = intensities.Sum();
            double weightedSumMz = 0;
            for (int i = 0; i < mzs.Count; i++)
            {
                weightedSumMz += mzs[i] * intensities[i];
            }
            X = weightedSumMz / Y;
        }

        #endregion Internal Methods

    }
}