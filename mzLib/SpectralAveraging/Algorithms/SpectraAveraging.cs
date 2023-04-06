using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MzLibUtil;

namespace SpectralAveraging;

public static class SpectraAveraging
{
    /// <summary>
    ///     Average a group of spectra in the jagged array format
    /// </summary>
    /// <param name="xArrays"> x values of spectra</param>
    /// <param name="yArrays">y values of spectra</param>
    /// <param name="parameters">Options for how to perform averaging</param>
    /// <returns>
    ///     Double jagged array where the first index is the x values and the second
    ///     is the y values for the averaged spectrum
    /// </returns>
    /// <exception cref="NotImplementedException">If merging type has not been implemented</exception>
    public static double[][] AverageSpectra(double[][] xArrays, double[][] yArrays,
        SpectralAveragingParameters parameters)
    {
        switch (parameters.SpectralAveragingType)
        {
            case SpectralAveragingType.MzBinning:
                return MzBinning(xArrays, yArrays, parameters);

            default:
                throw new MzLibException("Spectrum Averaging Type Not Yet Implemented");
        }
    }

    /// <summary>
    ///     Merges spectra into a two dimensional array of (m/z, int) values based upon their bin
    /// </summary>
    /// <param name="xArrays">xArrays of spectra to be averaged</param>
    /// <param name="yArrays">yArrays of spectra to be averaged</param>
    /// <param name="parameters">how to perform the averaging</param>
    /// <returns></returns>
    private static double[][] MzBinning(double[][] xArrays, double[][] yArrays,
        SpectralAveragingParameters parameters)
    {
        // get tics 
        var tics = yArrays.Select(p => p.Sum()).ToArray();
        var averageTic = tics.Average();

        // normalize spectra
        SpectraNormalization.NormalizeSpectra(yArrays, parameters.NormalizationType);

        // get bins
        var bins = GetBins(xArrays, yArrays, parameters.BinSize);

        // get weights
        var weights = SpectralWeighting.CalculateSpectraWeights(xArrays, yArrays, parameters.SpectralWeightingType);

        // reject outliers and average bins
        List<(double mz, double intensity)> averagedPeaks = new();
        Parallel.ForEach(Enumerable.Range(0, parameters.MaxThreadsToUsePerFile), (iterationIndex) =>
        {
            // each bin index that contains peaks
            var binIncidences = bins.Keys.ToList();

            // iterate through each bin index which contains peaks
            for (; iterationIndex < binIncidences.Count; iterationIndex += parameters.MaxThreadsToUsePerFile)
            {
                var peaksFromBin = bins[binIncidences[iterationIndex]];

                peaksFromBin = OutlierRejection.RejectOutliers(peaksFromBin, parameters);
                if (!peaksFromBin.Any()) continue;
                lock (averagedPeaks)
                    averagedPeaks.Add(AverageBin(peaksFromBin, weights));
            }
        });

        // return averaged
        var orderedAveragePeaks = averagedPeaks.Where(p => p.intensity != 0).OrderBy(p => p.mz).ToList();
        return new[]
        {
            orderedAveragePeaks.Select(p => p.mz).ToArray(),
            orderedAveragePeaks.Select(p =>
                parameters.NormalizationType == NormalizationType.AbsoluteToTic
                    ? p.intensity * averageTic
                    : p.intensity).ToArray()
        };
    }


    #region Helpers

    /// <summary>
    ///     Sorts spectra into bins
    /// </summary>
    /// <param name="xArrays">xArrays of spectra to be binned</param>
    /// <param name="yArrays">yArrays of spectra to be binned</param>
    /// <param name="binSize">size of bin</param>
    /// <returns>
    ///     Dictionary of bins where the key is the bin index and the
    ///     value is a list of all peaks in that particular bin
    /// </returns>
    /// <remarks>
    ///     There will be at minimum, one peak per bin per spectrum.
    ///     If a spectrum does not have a peak in a bin where all other spectra do,
    ///     then a peak with zero intensity will be added. This allows for efficient outlier rejection
    /// </remarks>
    private static Dictionary<int, List<BinnedPeak>> GetBins(double[][] xArrays, double[][] yArrays,
        double binSize)
    {
        var numSpectra = xArrays.Length;
        var minXValue = xArrays.MinBy(p => p.Min()).Min();

        // sort all spectra data into BinnedPeaks
        Dictionary<int, List<BinnedPeak>> bins = new();
        for (var i = 0; i < numSpectra; i++)
        for (var j = 0; j < xArrays[i].Length; j++)
        {
            var binIndex = (int)Math.Floor((xArrays[i][j] - minXValue) / binSize);
            var binValue = new BinnedPeak(binIndex, xArrays[i][j], yArrays[i][j], i);
            if (!bins.ContainsKey(binIndex))
                bins.Add(binIndex, new List<BinnedPeak> { binValue });
            else
                bins[binIndex].Add(binValue);
        }

        // add in zero values where a particular spectra did not have a peak in the bin
        foreach (var bin in bins)
        {
            var spectraInBin = bin.Value.Select(p => p.SpectraId).ToArray();
            var representativeBin = bin.Value.First();
            for (var i = 0; i < numSpectra; i++)
                if (!spectraInBin.Contains(i))
                    bin.Value.Add(new BinnedPeak(representativeBin.Bin, bin.Value.Average(p => p.Mz), 0, i));
        }

        return bins;
    }

    /// <summary>
    ///     Determines the weighted average of a list of binned peaks
    /// </summary>
    /// <param name="peaksInBin">peaks to average</param>
    /// <param name="weights">weights to average</param>
    /// <returns>tuple representing (average Mz value in bin, weighted average intensity in bin)</returns>
    private static (double, double) AverageBin(List<BinnedPeak> peaksInBin, Dictionary<int, double> weights)
    {
        double numerator = 0;
        double denominator = 0;

        foreach (var peak in peaksInBin)
        {
            numerator += peak.Intensity * weights[peak.SpectraId];
            denominator += weights[peak.SpectraId];
        }

        var mz = peaksInBin.Select(p => p.Mz).Average();
        var intensity = numerator / denominator;
        return (mz, intensity);
    }

    #endregion
}