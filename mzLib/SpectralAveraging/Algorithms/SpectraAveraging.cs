using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
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
    private static double[][] MzBinning(double[][] xArrays, double[][] yArrays, SpectralAveragingParameters parameters)
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
        var binIncidences = bins.Keys.ToList();
        (double mz, double intensity)[] averagedPeaks = new (double, double)[binIncidences.Count];
        var partitioner = Partitioner.Create(0, binIncidences.Count);

        Parallel.ForEach(partitioner, new ParallelOptions { MaxDegreeOfParallelism = parameters.MaxThreadsToUsePerFile }, (range, state) =>
        {
            var rejectionBuffers = new OutlierRejectionBuffers();
            for (int i = range.Item1; i < range.Item2; i++)
            {
                var peaksFromBin = bins[binIncidences[i]];

                int kept = OutlierRejection.RejectOutliers(CollectionsMarshal.AsSpan(peaksFromBin), parameters, rejectionBuffers);
                if (kept == 0) continue;
                peaksFromBin.RemoveRange(kept, peaksFromBin.Count - kept);

                averagedPeaks[i] = AverageBin(peaksFromBin, weights);
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
    internal static Dictionary<int, List<BinnedPeak>> GetBins(double[][] xArrays, double[][] yArrays,
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
            AddZeroIntensityPeaks(bin.Value, numSpectra);

        return bins;
    }

    /// <summary>
    ///     Adds a zero-intensity peak, at the bin's mean m/z, for each spectrum with no peak in the bin.
    /// </summary>
    /// <remarks>
    ///     The bin's own peaks are in spectrum order, so one pass finds the missing spectra. Each added peak's m/z is the
    ///     mean of everything in the bin at that point, the peaks added before it included; the m/z total is kept as the
    ///     bin grows rather than summed again for each added peak. It is the same sum - from the first peak to the last,
    ///     in order, as <see cref="Enumerable.Average{TSource}(IEnumerable{TSource}, Func{TSource, double})"/> sums -
    ///     divided by the same count, so each mean is the same double.
    /// </remarks>
    private static void AddZeroIntensityPeaks(List<BinnedPeak> bin, int numSpectra)
    {
        int ownPeaks = bin.Count;
        var representativeBin = bin[0];
        double mzTotal = bin[0].Mz;
        for (int p = 1; p < ownPeaks; p++)
            mzTotal += bin[p].Mz;

        int next = 0;
        for (var i = 0; i < numSpectra; i++)
        {
            while (next < ownPeaks && bin[next].SpectraId < i)
                next++;
            if (next < ownPeaks && bin[next].SpectraId == i)
                continue;

            double mz = mzTotal / bin.Count;
            bin.Add(new BinnedPeak(representativeBin.Bin, mz, 0, i));
            mzTotal += mz;
        }
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