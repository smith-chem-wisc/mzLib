using System.Collections.Generic;
using System.Linq;
using MzLibUtil;

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
                throw new MzLibException("Austin will fill this in when he returns from vacation");
            //return WeightByMrsNoiseEstimation(xArrays, yArrays);

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


    // NOTE: This will be implemented when Austin gets back from vacation
    // leaving him a nice place to put it

    /// <summary>
    ///     Given the noise estimates and the scale estimates, calculates the weight given to
    ///     each spectra when averaging using w_i = 1 / (k * noise_estimate)^2,
    ///     where k = scaleEstimate_reference / scaleEstimate_i
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
}