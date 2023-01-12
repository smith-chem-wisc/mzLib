using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;

namespace MzLibSpectralAveraging
{
    /// <summary>
    /// Weight each spectra 
    /// </summary>
    public static class SpectralWeighting
    {
        public static Dictionary<int, double> CalculateSpectraWeights(double[][] xArrays, double[][] yArrays, SpectraWeightingType spectraWeightingType)
        {
            switch (spectraWeightingType)
            {
                case SpectraWeightingType.WeightEvenly:
                    return WeightEvenly(xArrays.Length);

                case SpectraWeightingType.TicValue:
                    return WeightByTicValue(yArrays);

                case SpectraWeightingType.MrsNoiseEstimation:
                    return WeightByMrsNoiseEstimation(xArrays, yArrays);

                default:
                    throw new MzLibException("Spectra Weighting Type Not Implemented");
            }
        }


        private static Dictionary<int, double> WeightEvenly(int count)
        {
            var weights = new Dictionary<int, double>();
            for (int i = 0; i < count; i++)
            {
                weights.TryAdd(i, 1);
            }
            return weights;
        }

        private static Dictionary<int, double> WeightByTicValue(double[][] yArrays)
        {
            var weights = new Dictionary<int, double>();
            var tics = yArrays.Select(p => p.Sum()).ToArray();
            var maxTic = tics.Max();

            for (int i = 0; i < yArrays.Length; i++)
            {
                weights.TryAdd(i, tics[i] / maxTic);
            }
            return weights;
        }


        // NOTE: This will be implemented when Austin gets back from vacation

        /// <summary>
        /// Given the noise estimates and the scale estimates, calculates the weight given to
        /// each spectra when averaging using w_i = 1 / (k * noise_estimate)^2,
        /// where k = scaleEstimate_reference / scaleEstimate_i
        /// </summary>
        private static Dictionary<int, double> WeightByMrsNoiseEstimation(double[][] xArrays, double[][] yArrays)
        {
            var weights = new Dictionary<int, double>();
            // get noise and scale estimates
            //SortedDictionary<int, double> noiseEstimates = CalculateNoiseEstimates(pixelStack);
            //SortedDictionary<int, double> scaleEstimates = CalculateScaleEstimates(pixelStack);

            //// calculate weights
            //double referenceScale = scaleEstimates[0];
            //foreach (var entry in noiseEstimates)
            //{
            //    var successScale = scaleEstimates.TryGetValue(entry.Key,
            //        out double scale);
            //    if (!successScale) continue;

            //    var successNoise = noiseEstimates.TryGetValue(entry.Key,
            //        out double noise);
            //    if (!successNoise) continue;

            //    double k = referenceScale / scale;

            //    double weight = 1d / Math.Pow((k * noise), 2);

            //    weights.TryAdd(entry.Key, weight);
            //}

            return weights;
        }

        #region MrsNoiseHelpers

        ///// <summary>
        ///// Calculates noise estimates for each spectra using multi-resolution support noise estimation.
        ///// </summary>
        ///// <param name="waveletType">wavelet to be used in MRS noise estimation.</param>
        ///// <param name="epsilon">Noise estimate convergence to be reached before returning the noise estimate.</param>
        ///// <param name="maxIterations">Maximum number of iterations to be performed in the MRS noise estimation before returning.</param>
        //public static SortedDictionary<int, double> CalculateNoiseEstimates(PixelStack pixelStack, WaveletType waveletType = WaveletType.Haar,
        //    double epsilon = 0.01, int maxIterations = 25)
        //{
        //    ConcurrentDictionary<int, double> tempConcurrentDictionary = new();
        //    pixelStack.RecalculatedSpectra
        //        .Select((w, i) => new { Index = i, Array = w })
        //        .AsParallel()
        //        .ForAll(x =>
        //        {
        //            bool success = MRSNoiseEstimator.MRSNoiseEstimation(x.Array, epsilon, out double noiseEstimate,
        //                waveletType: waveletType, maxIterations: maxIterations);
        //            if (!success || double.IsNaN(noiseEstimate))
        //            {
        //                noiseEstimate = BasicStatistics.CalculateStandardDeviation(x.Array);
        //            }
        //            tempConcurrentDictionary.TryAdd(x.Index, noiseEstimate);
        //        });
        //    return new SortedDictionary<int, double>(tempConcurrentDictionary);
        //}

        ///// <summary>
        ///// Calculates the estimates of scale for each spectra in this object. Scale is determined by
        ///// taking the square root of the biweight midvariance. 
        ///// </summary>
        //public static SortedDictionary<int, double> CalculateScaleEstimates(PixelStack pixelStack)
        //{
        //    ConcurrentDictionary<int, double> tempScaleEstimates = new();
        //    pixelStack.RecalculatedSpectra
        //        .Select((w, i) => new { Index = i, Array = w })
        //        .AsParallel().ForAll(x =>
        //        {
        //            double scale = Math.Sqrt(BinweightMidvariance(x.Array));
        //            tempScaleEstimates.TryAdd(x.Index, Math.Sqrt(scale));
        //        });
        //    return new SortedDictionary<int, double>(tempScaleEstimates);
        //}

        ///// <summary>
        ///// Calcultes the biweight midvariance for an array. Algorithm orignally found here:
        ///// https://pixinsight.com/doc/tools/ImageIntegration/ImageIntegration.html#__equation_27__
        ///// </summary>
        ///// <param name="array">Array of doubles.</param>
        ///// <returns>The biweight midvariance.</returns>
        //private static double BinweightMidvariance(double[] array)
        //{
        //    double[] y_i = new double[array.Length];
        //    double[] a_i = new double[array.Length];
        //    double MAD_X = MedianAbsoluteDeviationFromMedian(array);
        //    double median = BasicStatistics.CalculateMedian(array);
        //    for (int i = 0; i < y_i.Length; i++)
        //    {
        //        y_i[i] = (array[i] - median) / (9d * MAD_X);
        //        if (y_i[i] < 1d)
        //        {
        //            a_i[i] = 1d;
        //        }
        //        else
        //        {
        //            a_i[i] = 0;
        //        }
        //    }

        //    // biweight midvariance calculation

        //    double denomSum = 0;
        //    double numeratorSum = 0;
        //    for (int i = 0; i < y_i.Length; i++)
        //    {
        //        numeratorSum += a_i[i] * Math.Pow(array[i] - median, 2) * Math.Pow(1 - y_i[i] * y_i[i], 4);
        //        denomSum += a_i[i] * (1 - 5 * y_i[i] * y_i[i]) * (1 - y_i[i] * y_i[i]);
        //    }

        //    return (double)y_i.Length * numeratorSum / Math.Pow(Math.Abs(denomSum), 2);
        //}

        ///// <summary>
        ///// Calculates the median absolute deviation from the median, which is then used in
        ///// calculating the biweight midvariance. Original algorithm found here:
        ///// https://pixinsight.com/doc/tools/ImageIntegration/ImageIntegration.html#__equation_26__
        ///// </summary>
        ///// <param name="array">Array of values to be calculated.</param>
        ///// <returns>The median absolute deviation from median.</returns>
        //private static double MedianAbsoluteDeviationFromMedian(double[] array)
        //{
        //    double arrayMedian = BasicStatistics.CalculateMedian(array);
        //    double[] results = new double[array.Length];
        //    for (int j = 0; j < array.Length; j++)
        //    {
        //        results[j] = Math.Abs(array[j] - arrayMedian);
        //    }

        //    return BasicStatistics.CalculateMedian(results);
        //}

        #endregion

    }

  
}
