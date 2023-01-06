using System;

namespace MzLibSpectralAveraging
{
    /// <summary>
    /// Class for processing signals by wavelet transforms. Currently,
    /// only the forward transform of the Maximum Overlap Discrete Wavelet Transform is implemented.
    /// </summary>
    /// <remarks>Note that the phrase "filter coefficients" has a distinct meaning separate from "wavelet coefficient" or "scaling coefficients".
    /// Filter coefficients are constants that refer to the wavelet or scaling values associated with a particular wavelet, e.g. a Haar wavelet.
    /// Wavelet coefficients or scaling coefficients refer to the output of a wavelet transform. </remarks>
    public static class WaveletMath
    {
        /// <summary>
        /// Calculates a forward maximum overlap discrete wavelet transform for a given signal at THE
        /// INDIVIDUAL SCALE. Outputs the wavelet coefficients and scaling coefficients for the calculated scale. Function was
        /// originally a C function written for an R package found at https://github.com/cran/wavelets/blob/master/src/modwt_forward.c. 
        /// </summary>
        /// <param name="V">Original signal</param>
        /// <param name="N">Length of original signal</param>
        /// <param name="j">Current scale</param>
        /// <param name="h">Wavelet filter array</param>
        /// <param name="g">Scaling filter array</param>
        /// <param name="L">Length of either wavelet or scaling filter (both are equivalent)</param>
        /// <param name="Wj">Resulting wavelet coefficients after performing the MODWT</param>
        /// <param name="Vj">Resulting Scaling coefficients after performing the MODWT</param>
        internal static void ModwtForward(double[] V, int N, int j,
            double[] h, double[] g, int L, ref double[] Wj, ref double[] Vj)
        {
            int t, k, n;
            double k_div;

            for (t = 0; t < N; t++)
            {
                k = t;
                Wj[t] = h[0] * V[k];
                Vj[t] = g[0] * V[k];
                for (n = 1; n < L; n++)
                {
                    k -= (int)Math.Pow(2, (j-1));
                    k_div = -(double)k / (double)N;
                    if (k < 0) k += (int)Math.Ceiling(k_div) * N;
                    Wj[t] += h[n] * V[k];
                    Vj[t] += g[n] * V[k]; 
                }
            }
        }

        /// <summary>
        /// The number of scales calculated in the MODWT depends on the length of the
        /// signal and the wavelet filter. This method calculates the number of scales to be calculated
        /// for the MODWT. 
        /// </summary>
        /// <param name="signalLength">Integer length of the input signal</param>
        /// <param name="filterLength">Integer length of the wavelet or scaling filter</param>
        /// <returns></returns>
        internal static int CalculateNumberScales(int signalLength, int filterLength)
        {
            double val = Math.Log(((signalLength - 1) / (filterLength - 1)) + 1) / Math.Log(2); 
            return (int)Math.Floor(val); 
        }

        /// <summary>
        /// Given a signal, wavelet filter, scaling filter and wavelet type,
        /// this method calculates the number of scales to be used in a MODWT transform
        /// and then calculates the MODWT at each scale.
        /// </summary>
        /// <remarks>The implementation of this MODWT uses a pyramidal algorithm. Each level calculation affects the next higher level.
        /// This is not the fastest possible implementation of a MODWT, but it is still very fast. If you want to make further optimizations,
        /// the correct route is to implement convolution by fourier transform and elementwise multiplication, and not parallelizing the for loop in this
        /// method.</remarks>
        /// <seealso cref="ModWtOutput"/>
        /// <seealso cref="ModwtForward"/>
        /// <param name="signal">Y-axis values. Typically intensity values. </param>
        /// <param name="waveletFilter">The wavelet filter coefficients.</param>
        /// <param name="scalingFilter">The scaling filter coefficients.</param>
        /// <param name="waveletType">The type of wavelet used in the transform.</param>
        internal static ModWtOutput ModWt(double[] signal, double[] waveletFilter, 
            double[] scalingFilter, WaveletType waveletType)
        {
            // calculate the number of scales to iterate over
            int numScales = CalculateNumberScales(signal.Length, waveletFilter.Length); 

            // use reflected boundary
            double[] reflectedSignal = CreateReflectedArray(signal); 

            var output = new ModWtOutput(numScales);

            double[] scalingCoeffs = new double[reflectedSignal.Length];
            reflectedSignal.CopyTo(scalingCoeffs, 0);


            for (int i = 0; i < numScales; i++)
            {
                double[] tempScalingCoeffs = new double[reflectedSignal.Length]; 
                scalingCoeffs.CopyTo(tempScalingCoeffs,0);

                double[] waveletCoeffs = new double[reflectedSignal.Length];
                ModwtForward(scalingCoeffs, reflectedSignal.Length, i + 1, waveletFilter,
                    scalingFilter, waveletFilter.Length, ref waveletCoeffs, 
                    ref tempScalingCoeffs);
                output.AddLevel(waveletCoeffs, scalingCoeffs, i, BoundaryType.Reflection, signal.Length, waveletFilter.Length);

                scalingCoeffs = tempScalingCoeffs; 
            }
            return output; 
        }
        /// <summary>
        /// Appends a reflection of the signal on the original array. 
        /// </summary>
        /// <param name="original"></param>
        /// <returns>A double array containing the original signal and a reflection of the original signal. Output size is
        /// 2x the original signal length. </returns>
        internal static double[] CreateReflectedArray(double[] original)
        {
            double[] reflectedArray = new double[original.Length*2];
            double[] copyOfOriginal = new double[original.Length];

            // copy original into new
            Buffer.BlockCopy(original, 0, copyOfOriginal, 0, sizeof(double)*copyOfOriginal.Length);
            // reverse copy of the original array
            Array.Reverse(copyOfOriginal);
            // Combine the original and the reverse arrays 
            Buffer.BlockCopy(original, 0, reflectedArray, 0, original.Length*sizeof(double));

            int reverseArrayOffset = sizeof(double) * (copyOfOriginal.Length); 
            // original array is copied starting at element zero
            Buffer.BlockCopy(copyOfOriginal, 0, 
                // dst array contains the original signal, so need to offset the copying 
                // of the reflected signal. 
                reflectedArray, reverseArrayOffset, 
                copyOfOriginal.Length * sizeof(double));
            return reflectedArray; 

        }
        /// <summary>
        /// Override method for the MODWT that takes the original signal and a WaveletFilter object. 
        /// </summary>
        /// <seealso cref="ModWt(double[],double[],double[],MzLibSpectralAveraging.WaveletType)"/>
        /// <seealso cref="WaveletFilter"/>
        /// <seealso cref="ModWtOutput"/>
        /// <param name="signal">The input signal. Typically an array of intensity values.</param>
        /// <param name="filter">WaveletFilter object that contains the wavelet and scaling filter coefficients.</param>
        /// <returns>A ModWtOutput object.</returns>
        internal static ModWtOutput ModWt(double[] signal, WaveletFilter filter)
        {
            return ModWt(signal, filter.WaveletCoefficients, filter.ScalingCoefficients, filter.WaveletType);
        }
    }
}
