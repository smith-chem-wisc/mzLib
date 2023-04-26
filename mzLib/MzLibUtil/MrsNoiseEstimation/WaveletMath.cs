using System;

namespace MzLibUtil.NoiseEstimation
{
    public class WaveletMath
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="V">Original signal</param>
        /// <param name="N">Length of original signal</param>
        /// <param name="j">Current scale</param>
        /// <param name="h">Wavelet filter</param>
        /// <param name="g">Scaling filter</param>
        /// <param name="L">Length of filter</param>
        /// <param name="Wj">Wavelet coefficients out</param>
        /// <param name="Vj">Scaling coefficients out</param>
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

        internal static int CalculateNumberScales(int signalLength, int filterLength)
        {
            double val = Math.Log(((signalLength - 1) / (filterLength - 1)) + 1) / Math.Log(2); 
            return (int)Math.Floor(val); 
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="signal"></param>
        /// <param name="waveletFilter"></param>
        /// <param name="scalingFilter"></param>
        /// <param name="numScales"></param>
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

        internal static ModWtOutput ModWt(double[] signal, WaveletFilter filters)
        {
            return ModWt(signal, filters.WaveletCoefficients, filters.ScalingCoefficients, filters.WaveletType);
        }
    }
}
