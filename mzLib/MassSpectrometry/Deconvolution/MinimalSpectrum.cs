using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry.Deconvolution
{
    // Consider defining this as a struct to increase performance
    public class MinimalSpectrum
    {
        public readonly double[] MzArray;
        public readonly double[] IntensityArray;
        public readonly double MostAbundantMz;
        public readonly int Charge;

        public MinimalSpectrum(double[] mzArray, double[] intensityArray, int charge = 0)
        {
            MzArray = mzArray;
            IntensityArray = intensityArray;
            MostAbundantMz = GetMostAbundantMz(mzArray, intensityArray);
            Charge = charge;
        }

        internal double[] GetMzs()
        {
            double[] mzArrayCopy = new double[MzArray.Length];
            Array.Copy(MzArray, mzArrayCopy, MzArray.Length);
            return mzArrayCopy;
        }

        internal double[] GetIntensities()
        {
            double[] intensityArrayCopy = new double[IntensityArray.Length];
            Array.Copy(MzArray, intensityArrayCopy, IntensityArray.Length);
            return intensityArrayCopy;
        }

        /// <summary>
        /// Returns the charge, or 0 if charge was not assigned
        /// </summary>
        /// <returns></returns>
        internal int GetCharge()
        {
            return Charge;
        }

        internal static double GetMostAbundantMz(double[] mzArray, double[] intensityArray)
        {
            double mostAbundantMz = 0;
            double maxIntensity = 0;
            for (int i = 0; i < mzArray.Length; i++)
            {
                if (intensityArray[i] > maxIntensity)
                {
                    maxIntensity = intensityArray[i];
                    mostAbundantMz = mzArray[i];
                }
            }

            return mostAbundantMz;

        }
    }
}
