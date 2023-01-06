using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MzLibSpectralAveraging
{
    public static class SpectraWeighting
    {
        public static SortedDictionary<int, double> CalculateSpectraWeights(BinnedSpectra binnedSpectra, SpectraWeightingType weightingType)
        {
            switch (weightingType)
            {
                case SpectraWeightingType.None:
                    return WeightEvenly(binnedSpectra);

                case SpectraWeightingType.TicValue:
                    return WeightByTicValue(binnedSpectra);

                case SpectraWeightingType.MrsNoiseEstimation:
                    return WeightByMrsNoiseEstimation(binnedSpectra);

                default:
                    throw new NotImplementedException("Spectra Weighting Type Not Implemented");
            }
        }

        private static SortedDictionary<int, double> WeightEvenly(BinnedSpectra binnedSpectra)
        {
            var weights = new SortedDictionary<int, double>();

            for (int i = 0; i < binnedSpectra.NumSpectra; i++)
            {
                weights.TryAdd(i, 1);
            }
            return weights;
        }

        private static SortedDictionary<int, double> WeightByTicValue(BinnedSpectra binnedSpectra)
        {
            var weights = new SortedDictionary<int, double>();
            var maxTic = binnedSpectra.Tics.Max();

            for (int i = 0; i < binnedSpectra.NumSpectra; i++)
            {
                weights.TryAdd(i, binnedSpectra.Tics[i] / maxTic);
            }
            return weights;
        }

        private static SortedDictionary<int, double> WeightByMrsNoiseEstimation(BinnedSpectra binnedSpectra)
        {
            var weights = new SortedDictionary<int, double>();
            // get noise and scale estimates
            SortedDictionary<int, double> noiseEstimates = CalculateNoiseEstimates(binnedSpectra);
            SortedDictionary<int, double> scaleEstimates = CalculateScaleEstimates(binnedSpectra);

            // calculate weights
            double referenceScale = scaleEstimates[binnedSpectra.ReferenceSpectra];
            foreach (var entry in noiseEstimates)
            {
                var successScale = scaleEstimates.TryGetValue(entry.Key,
                    out double scale);
                if (!successScale) continue;

                var successNoise = noiseEstimates.TryGetValue(entry.Key,
                    out double noise);
                if (!successNoise) continue;

                double k = referenceScale / scale;

                double weight = 1d / Math.Pow((k * noise), 2);

                weights.TryAdd(entry.Key, weight);
            }

            return weights;
        }

        #region MrsNoiseHelpers

        /// <summary>
        /// Calculates noise estimates for each spectra using multi-resolution support noise estimation.
        /// </summary>
        /// <param name="waveletType">wavelet to be used in MRS noise estimation.</param>
        /// <param name="epsilon">Noise estimate convergence to be reached before returning the noise estimate.</param>
        /// <param name="maxIterations">Maximum number of iterations to be performed in the MRS noise estimation before returning.</param>
        public static SortedDictionary<int, double> CalculateNoiseEstimates(BinnedSpectra binnedSpectra, WaveletType waveletType = WaveletType.Haar,
            double epsilon = 0.01, int maxIterations = 25)
        {
            ConcurrentDictionary<int, double> tempConcurrentDictionary = new();
            binnedSpectra.RecalculatedSpectra
                .Select((w, i) => new { Index = i, Array = w })
                .AsParallel()
                .ForAll(x =>
                {
                    bool success = MRSNoiseEstimator.MRSNoiseEstimation(x.Array, epsilon, out double noiseEstimate,
                        waveletType: waveletType, maxIterations: maxIterations);
                    if (!success || double.IsNaN(noiseEstimate))
                    {
                        noiseEstimate = BasicStatistics.CalculateStandardDeviation(x.Array);
                    }
                    tempConcurrentDictionary.TryAdd(x.Index, noiseEstimate);
                });
            return new SortedDictionary<int, double>(tempConcurrentDictionary);
        }

        /// <summary>
        /// Calculates the estimates of scale for each spectra in this object. Scale is determined by
        /// taking the square root of the biweight midvariance. 
        /// </summary>
        public static SortedDictionary<int, double> CalculateScaleEstimates(BinnedSpectra binnedSpectra)
        {
            ConcurrentDictionary<int, double> tempScaleEstimates = new();
            binnedSpectra.RecalculatedSpectra
                .Select((w, i) => new { Index = i, Array = w })
                .AsParallel().ForAll(x =>
                {
                    double scale = Math.Sqrt(BinweightMidvariance(x.Array));
                    tempScaleEstimates.TryAdd(x.Index, Math.Sqrt(scale));
                });
            return new SortedDictionary<int, double>(tempScaleEstimates);
        }

        /// <summary>
        /// Calcultes the biweight midvariance for an array. Algorithm orignally found here:
        /// https://pixinsight.com/doc/tools/ImageIntegration/ImageIntegration.html#__equation_27__
        /// </summary>
        /// <param name="array">Array of doubles.</param>
        /// <returns>The biweight midvariance.</returns>
        private static double BinweightMidvariance(double[] array)
        {
            double[] y_i = new double[array.Length];
            double[] a_i = new double[array.Length];
            double MAD_X = MedianAbsoluteDeviationFromMedian(array);
            double median = BasicStatistics.CalculateMedian(array);
            for (int i = 0; i < y_i.Length; i++)
            {
                y_i[i] = (array[i] - median) / (9d * MAD_X);
                if (y_i[i] < 1d)
                {
                    a_i[i] = 1d;
                }
                else
                {
                    a_i[i] = 0;
                }
            }

            // biweight midvariance calculation

            double denomSum = 0;
            double numeratorSum = 0;
            for (int i = 0; i < y_i.Length; i++)
            {
                numeratorSum += a_i[i] * Math.Pow(array[i] - median, 2) * Math.Pow(1 - y_i[i] * y_i[i], 4);
                denomSum += a_i[i] * (1 - 5 * y_i[i] * y_i[i]) * (1 - y_i[i] * y_i[i]);
            }

            return (double)y_i.Length * numeratorSum / Math.Pow(Math.Abs(denomSum), 2);
        }

        /// <summary>
        /// Calculates the median absolute deviation from the median, which is then used in
        /// calculating the biweight midvariance. Original algorithm found here:
        /// https://pixinsight.com/doc/tools/ImageIntegration/ImageIntegration.html#__equation_26__
        /// </summary>
        /// <param name="array">Array of values to be calculated.</param>
        /// <returns>The median absolute deviation from median.</returns>
        private static double MedianAbsoluteDeviationFromMedian(double[] array)
        {
            double arrayMedian = BasicStatistics.CalculateMedian(array);
            double[] results = new double[array.Length];
            for (int j = 0; j < array.Length; j++)
            {
                results[j] = Math.Abs(array[j] - arrayMedian);
            }

            return BasicStatistics.CalculateMedian(results);
        }

        #endregion

    }
}
