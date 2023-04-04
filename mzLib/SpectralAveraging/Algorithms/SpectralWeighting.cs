using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using MzLibUtil.NoiseEstimation;

namespace SpectralAveraging;

/// <summary>
///     Weight each spectra
/// </summary>
public static class SpectralWeighting
{
    /// <summary>
    ///     Calls the specific weighting function to determine weight to be applied for each spectra
    /// </summary>
    /// <param name="xArrays">xArrays of spectra to determine weights from</param>
    /// <param name="yArrays">yArrays of spectra to determine weights from</param>
    /// <param name="spectraWeightingType">how to weight the spectra</param>
    /// <returns>Dictionary of weights where the key is the spectra index and the value is the weight</returns>
    /// <exception cref="MzLibException"></exception>
    public static Dictionary<int, double> CalculateSpectraWeights(double[][] xArrays, double[][] yArrays,
        SpectraWeightingType spectraWeightingType)
    {
        switch (spectraWeightingType)
        {
            case SpectraWeightingType.WeightEvenly:
                return WeightEvenly(xArrays.Length);

            case SpectraWeightingType.TicValue:
                return WeightByTicValue(yArrays);

            case SpectraWeightingType.MrsNoiseEstimation:
                return WeightByMrsNoiseEstimation(yArrays);

            default:
                throw new MzLibException("Spectra Weighting Type Not Implemented");
        }
    }

    /// <summary>
    ///     Weights each spectra evenly
    /// </summary>
    /// <param name="count">number of spectra</param>
    /// <returns></returns>
    private static Dictionary<int, double> WeightEvenly(int count)
    {
        var weights = new Dictionary<int, double>();
        for (var i = 0; i < count; i++) weights.TryAdd(i, 1);
        return weights;
    }

    /// <summary>
    ///     Weight relative to the maximum tic value
    /// </summary>
    /// <param name="yArrays"></param>
    /// <returns></returns>
    private static Dictionary<int, double> WeightByTicValue(double[][] yArrays)
    {
        var weights = new Dictionary<int, double>();
        var tics = yArrays.Select(p => p.Sum()).ToArray();
        var maxTic = tics.Max();

        for (var i = 0; i < yArrays.Length; i++) weights.TryAdd(i, tics[i] / maxTic);
        return weights;
    }

    /// <summary>
    /// Given the noise estimates and the scale estimates, calculates the weight given to
    /// each spectra when averaging using w_i = 1 / (k * noise_estimate)^2,
    /// where k = scaleEstimate_reference / scaleEstimate_i. Reference spectra defaults to the
    /// first spectra
    /// </summary>
    private static Dictionary<int, double> WeightByMrsNoiseEstimation(double[][] yArrays)
    {
        var weights = new Dictionary<int, double>();
        // get noise estimate
        ConcurrentDictionary<int, double> noiseEstimates = new(); 

        yArrays
            .Select((w,i) => new {Index = i, Array = w})
            .AsParallel()
            .ForAll(x =>
            {
                bool mrsSuccess = MRSNoiseEstimator.MRSNoiseEstimation(x.Array, 0.01, out double noiseEstimate);
                if (!mrsSuccess || double.IsNaN(noiseEstimate) || noiseEstimate == 0d)
                {
                    noiseEstimate = x.Array.StandardDeviation();
                }

                noiseEstimates.TryAdd(x.Index, noiseEstimate);
            });
        // get scale estimate 
        ConcurrentDictionary<int, double> scaleEstimates = new(); 
        yArrays.Select((w,i) => new {Index = i, Array = w})
            .AsParallel()
            .ForAll(x =>
            {
                double scale = Math.Sqrt(BasicStatistics.BiweightMidvariance(x.Array));
                scaleEstimates.TryAdd(x.Index, Math.Sqrt(scale));
            });

        CalculateWeights(noiseEstimates, scaleEstimates, weights);

        return weights;
    }
    /// <summary>
    /// Given the noise estimates and the scale estimates, calculates the weight given to
    /// each spectra when averaging using w_i = 1 / (k * noise_estimate)^2,
    /// where k = scaleEstimate_reference / scaleEstimate_i
    /// </summary>
    private static void CalculateWeights(IDictionary<int, double> noiseEstimates, 
        IDictionary<int, double> scaleEstimates, IDictionary<int, double> weights,int referenceSpectra = 0)
    {
        double referenceScale = scaleEstimates[referenceSpectra];
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
    }

    
}