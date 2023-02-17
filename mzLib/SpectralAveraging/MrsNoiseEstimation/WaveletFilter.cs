using System;

namespace SpectralAveraging
{
    public class WaveletFilter
    {
        public double[] WaveletCoefficients { get; private set; }
        public double[] ScalingCoefficients { get; private set; }
        public WaveletType WaveletType { get; private set; }

        private void CreateFiltersFromCoeffs(double[] filterCoeffs)
        {
            WaveletCoefficients = new double[filterCoeffs.Length];
            ScalingCoefficients = new double[filterCoeffs.Length];

            // calculate wavelet coefficients
            for (int i = 0; i < ScalingCoefficients.Length; i++)
            {
                ScalingCoefficients[i] = filterCoeffs[i] / Math.Sqrt(2d);
            }
            WaveletCoefficients = WaveletMathUtils.QMF(ScalingCoefficients, inverse: true);
        }

        public void CreateFiltersFromCoeffs(WaveletType waveletType)
        {
            switch (waveletType)
            {
                case WaveletType.Haar:
                {
                    WaveletType = WaveletType.Haar;
                    CreateFiltersFromCoeffs(_haarCoefficients);
                    return;
                }
                //case WaveletType.Db4:
                //{
                //    WaveletType = WaveletType.Db4; 
                //    CreateFiltersFromCoeffs(_db4Coefficients);
                //    return; 
                //}
            }
        }
        #region Wavelet Coefficients
        private readonly double[] _haarCoefficients =
        {
            0.7071067811865475,
            0.7071067811865475
        };

        private readonly double[] _db4Coefficients =
        {
            0.2304, 
            0.7148, 
            0.6309, 
            -0.0280, 
            -0.1870, 
            0.0308, 
            0.0329, 
            -0.0106
        };
        #endregion
    }
}