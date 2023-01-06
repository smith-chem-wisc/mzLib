using System;
using System.Collections.Generic;
using System.Linq;

namespace MzLibSpectralAveraging
{
    /// <summary>
    /// Class used to calculate Multi Resolution Support noise estimation. Algorithm is based on this paper:
    /// Jean-Luc Starck and Fionn Murtagh (1998), Automatic Noise Estimation from the Multiresolution Support, Publications of the Royal Astronomical Society of the Pacific, vol. 110, pp. 193–199. 
    /// </summary>
    public static class MRSNoiseEstimator
    {
        /// <summary>
        /// Method to calculate the estimate of the noise standard deviation using the Multi-Resolution
        /// Support methods defined by Jean‐Luc Starck and Fionn Murtagh (Feb. 1998). Calculates a Maximum Overlap
        /// Discrete wavelet transform, then creates a multi-resolution support to find a set of noise pixels. Takes the standard
        /// deviation of these only the noise pixels to use in estimating the standard deviation of the noise. 
        /// </summary>
        /// <seealso cref="BooleanizeLevel"/>
        /// <seealso cref="CreateMultiResolutionSupport"/>
        /// <remarks>The Haar wavelet is the default, and I don't see a lot of reason that it would need to be different. It will be the fastest wavelet to use, as it is only two points.
        /// Further, the other wavelet that I've implemented, the Db4, has similar estimates in noise, but is significantly slower. Since we're dealing with spectra that will include many points,
        /// I think the speed is an important advantage supporting the Haar wavelet as the default.</remarks>
        /// <param name="signal">A time or frequency domain signal.</param>
        /// <param name="epsilon">Determines threshold for convergence between iterations.</param>
        /// <param name="noiseEstimate">The noise estimate value for a the given signal.</param>
        /// <param name="maxIterations">Maximum number of iterations to perform.</param>
        /// <param name="waveletType">Wavelet used in creating the multi-resolution support. Default is Haar, and I don't see any reason to make it different.</param>
        /// <returns>A boolean value indicating if the algorithm converged on a value or not. Used to determine
        /// if alternate noise estimation is required.</returns>
        public static bool MRSNoiseEstimation(double[] signal, double epsilon,
            out double noiseEstimate, int maxIterations = 25, 
            WaveletType waveletType = WaveletType.Haar)
        {
            int iterations = 0; 
            // 1. Estimate the standard deviation of the noise in the original signal. 
            double stdevPrevious = BasicStatistics.CalculateStandardDeviation(signal);

            // 2. Compute the modwt of the image
            WaveletFilter filter = new(); 
            filter.CreateFiltersFromCoeffs(waveletType);
            ModWtOutput wtOutput = WaveletMath.ModWt(signal, filter);

            // 3. Set n = 0; (this is done implicitly while in the do while loop.)

            double[] signalIterable = new double[signal.Length]; 
            signal.CopyTo(signalIterable, 0);


            signalIterable = CreateSmoothedSignal(signalIterable, wtOutput);
            
            double criticalVal;
            double stdevNext;
            do
            {
                // 4. Compute the multiresolution support M that is derived from the wavelet coefficients
                // and the standard deviation of the noise at each level. 
                // 5. Select all points that belong to the noise; they don't have an significant coefficients above noise 
                var booleanizedLevels = BooleanizeLevels(wtOutput,
                    stdevPrevious, 1.97); 
                int[] mrsIndices = CreateMultiResolutionSupport(booleanizedLevels);

                // 6. For the selected pixels, calculate original array - smoothed array and compute the standard deviation 
                // for those values. 
                // don't modify the original signal, use a deep copy instead: 
                stdevNext = ComputeStdevOfNoisePixels(signalIterable, mrsIndices);

                // 7. n = n + l. 
                // 8. start again at 4 if sigma_I^n - sigma_I^(n-1) / sigma_I^(n) > epsilon. 
                criticalVal = Math.Abs(stdevNext - stdevPrevious) / stdevPrevious;

                // setup for next iteration 
                iterations++; 
                stdevPrevious = stdevNext;
            } while (criticalVal > epsilon && iterations <= maxIterations);

            noiseEstimate = stdevPrevious; 
            return iterations < maxIterations;
        }

        /// <summary>
        /// Converts the ModWtOutput into a list of int arrays based on a noise estimate. Called once per iteration to determine
        /// what is a noise pixel versus a signal pixel. 
        /// </summary>
        /// <param name="wtOutput">Output object of the Mod wavelet transform.</param>
        /// <param name="noiseEstimate">The initial noise estimation of the signal given for this iteration.</param>
        /// <param name="noiseThreshold">A multiple of the noise estimate. noiseThreshold * noiseEstimate
        /// is the cutoff value determining whether or not a pixel is signal or noise.</param>
        /// <returns>List<int[]>. Each element in the int[] is either a 0 or positive integer, representing a noise value (if 0)
        /// or a signal value (if int[i] > 0).</int></returns>
        private static List<int[]> BooleanizeLevels(ModWtOutput wtOutput, double noiseEstimate, double noiseThreshold)
        {
            List<int[]> booleanizedLevels = new();
            foreach (Level level in wtOutput.Levels)
            {
                booleanizedLevels.Add(BooleanizeLevel(level, noiseEstimate, noiseThreshold));
            }
            return booleanizedLevels;
        }

        /// <summary>
        /// Private method used to booleanize individual levels. For each pixel in the Level to be booleanized,
        /// determines if the values of that pixel is singificant, i.e. > threshold. If the pixel is significant,
        /// a one is added at that position to the output int[] array. 
        /// </summary>
        /// <param name="level">The level within a ModWtOuput.</param>
        /// <param name="noiseEstimate">Initial noise estimate of the signal.</param>
        /// <param name="threshold">Multiple of the noise estimate used to define the noise cutoff level.</param>
        /// <returns>An integer array where each index corresponds to a position in the signal. If the signal at that position is significant
        /// at the given level, the value at the index is >0.</returns>
        private static int[] BooleanizeLevel(Level level, double noiseEstimate, double threshold)
        {
            int[] results = new int[level.WaveletCoeff.Length];
            double valToExceed = threshold * noiseEstimate;

            for (int i = 0; i < results.Length; i++)
            {
                if (Math.Abs(level.WaveletCoeff[i]) >= valToExceed)
                {
                    results[i] = 1;
                }
                else
                {
                    results[i] = 0;
                }
            }
            return results;
        }
        /// <summary>
        /// Returns a signal that is smoothed by subtracting the summed wavelet coefficients from the original signal. Used in each iteration of MRS noise estimation process.
        /// Signal should be essentially noiseless. 
        /// </summary>
        /// <param name="originalSignal">The original signal before modwt.</param>
        /// <param name="output">The ModWtOutput object containing the wavelet coefficients from the modwt of the original signal.</param>
        /// <returns>The original signal - the sum of the wavelet coefficients</returns>
        private static double[] CreateSmoothedSignal(double[] originalSignal, ModWtOutput output)
        {
            double[] results = new double[originalSignal.Length];
            int lastLevel = output.MaxScale; 
            double[] summedWt = output.Levels[lastLevel-1].ScalingCoeff;

            for (int i = 0; i < results.Length; i++)
            {
                results[i] = originalSignal[i] - summedWt[i]; 
            }

            return results; 
        }

        /// <summary>
        /// Sums the booleanized levels (a list of int[]) to creates a single integer array. 
        /// </summary>
        /// <param name="booleanizedLevels"></param>
        /// <returns>An integer array where 0 indicates a noise pixel and a positive integer >0 equals a
        /// signal pixel.</returns>
        private static int[] CreateMultiResolutionSupport(List<int[]> booleanizedLevels)
        {
            int[] outputArray = new int[booleanizedLevels[0].Length];
            for (int i = 0; i < outputArray.Length; i++)
            {
                outputArray[i] = booleanizedLevels.Select(k => k.ElementAt(i)).Sum();
            }
            return outputArray;
        }

        /// <summary>
        /// Computes the standard deviation of the noise pixels of a signal given the signal, the
        /// output of the ModWt wavelet transform and the indexes generated from multi-resolution support. 
        /// </summary>
        /// <param name="signal">A signal that has been smoothed by subtracting the wavelet coefficients from the original signal.</param>
        /// <param name="noiseIndices">An array of noise indices with either a zero or a positive integer value.
        /// Zero represents a noise pixel; positive integer represents a non-noise pixel.</param>
        /// <returns>The standard deviation of the noise pixels.</returns>
        private static double ComputeStdevOfNoisePixels(double[] signal, int[] noiseIndices)
        {
            List<double> noiseValues = new();
            // for each index in noiseIndices, get the values of the noise index, which will be either zero or one. 
            // If the value is zero, that means the value is a noise value. Otherwise, the index represents a signal. 
            for (int k = 0; k < noiseIndices.Length; k++)
            {
                if (noiseIndices[k] == 0)
                {
                    noiseValues.Add(signal.ElementAt(k));
                }
            }
            // returns the standard deviation of the noise values. 
            return BasicStatistics.CalculateStandardDeviation(noiseValues);
        }
    }

    
}
