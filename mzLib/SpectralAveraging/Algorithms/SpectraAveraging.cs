using System;
using System.Buffers;
using System.Collections.Concurrent;
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
    private static double[][] MzBinning(double[][] xArrays, double[][] yArrays, SpectralAveragingParameters parameters)
    {
        // get tics 
        var tics = yArrays.Select(p => p.Sum()).ToArray();
        var averageTic = tics.Average();

        // normalize spectra
        SpectraNormalization.NormalizeSpectra(yArrays, parameters.NormalizationType);

        // get bins
        using var bins = GetBins(xArrays, yArrays, parameters.BinSize);

        // get weights
        var weights = new SpectrumWeights(
            SpectralWeighting.CalculateSpectraWeights(xArrays, yArrays, parameters.SpectralWeightingType), xArrays.Length);

        // reject outliers and average bins
        var averagedPeaks = ArrayPool<(double mz, double intensity)>.Shared.Rent(bins.Count);
        var partitioner = Partitioner.Create(0, bins.Count);

        Parallel.ForEach(partitioner, new ParallelOptions { MaxDegreeOfParallelism = parameters.MaxThreadsToUsePerFile }, (range, state) =>
        {
            var rejectionBuffers = new OutlierRejectionBuffers();
            for (int i = range.Item1; i < range.Item2; i++)
            {
                var peaksFromBin = bins[i];

                int kept = OutlierRejection.RejectOutliers(peaksFromBin, parameters, rejectionBuffers);
                averagedPeaks[i] = kept == 0 ? default : AverageBin(peaksFromBin[..kept], weights);
            }
        });

        // return averaged
        var orderedAveragePeaks = new ArraySegment<(double mz, double intensity)>(averagedPeaks, 0, bins.Count)
            .Where(p => p.intensity != 0).OrderBy(p => p.mz).ToList();
        ArrayPool<(double mz, double intensity)>.Shared.Return(averagedPeaks);
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
    ///     The bins, numbered in the order they are first met, each holding all the peaks that fall in it
    /// </returns>
    /// <remarks>
    ///     There will be at minimum, one peak per bin per spectrum.
    ///     If a spectrum does not have a peak in a bin where all other spectra do,
    ///     then a peak with zero intensity will be added. This allows for efficient outlier rejection
    ///     <para>
    ///     All bins share one rented array (see <see cref="BinnedSpectra"/>) instead of a list each: a top-down MS1
    ///     window has tens of thousands of bins, and averaging a file bins every window. The bins come out in the order
    ///     a <c>Dictionary&lt;int, List&lt;BinnedPeak&gt;&gt;</c> keyed on the bin index enumerated them, holding the
    ///     same peaks in the same order.
    ///     </para>
    /// </remarks>
    internal static BinnedSpectra GetBins(double[][] xArrays, double[][] yArrays, double binSize)
    {
        var numSpectra = xArrays.Length;
        var minXValue = xArrays.MinBy(p => p.Min()).Min();

        int totalPeaks = 0;
        for (var i = 0; i < numSpectra; i++)
            totalPeaks += xArrays[i].Length;

        // the bin index of every peak
        var intPool = ArrayPool<int>.Shared;
        int[] binOfPeak = intPool.Rent(totalPeaks);
        int lowestIndex = int.MaxValue;
        int highestIndex = int.MinValue;
        int peak = 0;
        for (var i = 0; i < numSpectra; i++)
        for (var j = 0; j < xArrays[i].Length; j++)
        {
            var binIndex = (int)Math.Floor((xArrays[i][j] - minXValue) / binSize);
            binOfPeak[peak++] = binIndex;
            lowestIndex = Math.Min(lowestIndex, binIndex);
            highestIndex = Math.Max(highestIndex, binIndex);
        }

        // Number the bins in the order they are first met, and count each one's peaks and spectra; binOfPeak then
        // holds each peak's bin number. Bin indices are looked up in a table when they span a range not much wider
        // than the number of peaks, as they do for spectra binned near their own resolution, and in a dictionary
        // otherwise.
        int[] binIndexOfBin = intPool.Rent(totalPeaks);
        int[] ownPeaks = intPool.Rent(totalPeaks);
        int[] spectraInBin = intPool.Rent(totalPeaks);
        int[] lastSpectrum = intPool.Rent(totalPeaks);
        long indexRange = (long)highestIndex - lowestIndex + 1;
        int[] binOfIndex = null;
        Dictionary<int, int> binOfIndexLookup = null;
        if (indexRange <= 4L * totalPeaks + 1024)
        {
            binOfIndex = intPool.Rent((int)indexRange);
            Array.Fill(binOfIndex, -1, 0, (int)indexRange);
        }
        else
        {
            binOfIndexLookup = new Dictionary<int, int>();
        }

        int binCount = 0;
        peak = 0;
        for (var i = 0; i < numSpectra; i++)
        for (var j = 0; j < xArrays[i].Length; j++)
        {
            int binIndex = binOfPeak[peak];
            int bin;
            bool seen;
            if (binOfIndex is not null)
            {
                bin = binOfIndex[binIndex - lowestIndex];
                seen = bin >= 0;
                if (!seen)
                    binOfIndex[binIndex - lowestIndex] = bin = binCount;
            }
            else
            {
                seen = binOfIndexLookup.TryGetValue(binIndex, out bin);
                if (!seen)
                    binOfIndexLookup.Add(binIndex, bin = binCount);
            }

            if (!seen)
            {
                binCount++;
                binIndexOfBin[bin] = binIndex;
                ownPeaks[bin] = 0;
                spectraInBin[bin] = 0;
                lastSpectrum[bin] = -1;
            }

            ownPeaks[bin]++;
            if (lastSpectrum[bin] != i)
            {
                lastSpectrum[bin] = i;
                spectraInBin[bin]++;
            }

            binOfPeak[peak++] = bin;
        }

        // lay the bins out one after another, leaving room for a zero-intensity peak per spectrum missing from each
        int[] start = intPool.Rent(binCount);
        int[] length = intPool.Rent(binCount);
        int slots = 0;
        for (int bin = 0; bin < binCount; bin++)
        {
            start[bin] = slots;
            length[bin] = ownPeaks[bin] + numSpectra - spectraInBin[bin];
            slots += length[bin];
        }

        // sort all spectra data into BinnedPeaks, each bin's own peaks in spectrum order at its front
        BinnedPeak[] peaks = ArrayPool<BinnedPeak>.Shared.Rent(slots);
        int[] placed = lastSpectrum;
        Array.Clear(placed, 0, binCount);
        peak = 0;
        for (var i = 0; i < numSpectra; i++)
        for (var j = 0; j < xArrays[i].Length; j++)
        {
            int bin = binOfPeak[peak++];
            peaks[start[bin] + placed[bin]++] = new BinnedPeak(binIndexOfBin[bin], xArrays[i][j], yArrays[i][j], i);
        }

        // add in zero values where a particular spectra did not have a peak in the bin
        for (int bin = 0; bin < binCount; bin++)
            AddZeroIntensityPeaks(peaks.AsSpan(start[bin], length[bin]), ownPeaks[bin], numSpectra);

        intPool.Return(binOfPeak);
        intPool.Return(binIndexOfBin);
        intPool.Return(ownPeaks);
        intPool.Return(spectraInBin);
        intPool.Return(lastSpectrum);
        if (binOfIndex is not null)
            intPool.Return(binOfIndex);

        return new BinnedSpectra(peaks, start, length, binCount);
    }

    /// <summary>
    ///     Adds a zero-intensity peak, at the bin's mean m/z, for each spectrum with no peak in the bin. The first
    ///     <paramref name="ownPeaks"/> entries of <paramref name="bin"/> are the bin's own peaks; the rest are filled.
    /// </summary>
    /// <remarks>
    ///     The bin's own peaks are in spectrum order, so one pass finds the missing spectra. Each added peak's m/z is the
    ///     mean of everything in the bin at that point, the peaks added before it included; the m/z total is kept as the
    ///     bin grows rather than summed again for each added peak. It is the same sum - from the first peak to the last,
    ///     in order, as <see cref="Enumerable.Average{TSource}(IEnumerable{TSource}, Func{TSource, double})"/> sums -
    ///     divided by the same count, so each mean is the same double.
    /// </remarks>
    private static void AddZeroIntensityPeaks(Span<BinnedPeak> bin, int ownPeaks, int numSpectra)
    {
        int count = ownPeaks;
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

            double mz = mzTotal / count;
            bin[count++] = new BinnedPeak(representativeBin.Bin, mz, 0, i);
            mzTotal += mz;
        }
    }

    /// <summary>
    ///     Determines the weighted average of a list of binned peaks
    /// </summary>
    /// <param name="peaksInBin">peaks to average</param>
    /// <param name="weights">weights to average</param>
    /// <returns>tuple representing (average Mz value in bin, weighted average intensity in bin)</returns>
    /// <remarks>
    ///     The mean m/z is summed from the first peak to the last and divided by the count, as
    ///     <see cref="Enumerable.Average(IEnumerable{double})"/> computes it.
    /// </remarks>
    private static (double, double) AverageBin(ReadOnlySpan<BinnedPeak> peaksInBin, SpectrumWeights weights)
    {
        double numerator = 0;
        double denominator = 0;

        foreach (var peak in peaksInBin)
        {
            numerator += peak.Intensity * weights[peak.SpectraId];
            denominator += weights[peak.SpectraId];
        }

        double mzTotal = peaksInBin[0].Mz;
        for (int i = 1; i < peaksInBin.Length; i++)
            mzTotal += peaksInBin[i].Mz;

        var mz = mzTotal / peaksInBin.Length;
        var intensity = numerator / denominator;
        return (mz, intensity);
    }

    /// <summary>
    ///     The weight of each spectrum, read by its index. Every weighting method weighs every spectrum, so the weights
    ///     are copied into an array once rather than hashed for each peak of each bin; a spectrum without one is looked
    ///     up in the dictionary, which throws exactly as looking it up for each peak did.
    /// </summary>
    private sealed class SpectrumWeights
    {
        private readonly Dictionary<int, double> _weights;
        private readonly double[] _weightOf;
        private readonly bool[] _hasWeight;

        internal SpectrumWeights(Dictionary<int, double> weights, int numSpectra)
        {
            _weights = weights;
            _weightOf = new double[numSpectra];
            _hasWeight = new bool[numSpectra];
            for (int i = 0; i < numSpectra; i++)
                _hasWeight[i] = weights.TryGetValue(i, out _weightOf[i]);
        }

        internal double this[int spectrum] => _hasWeight[spectrum] ? _weightOf[spectrum] : _weights[spectrum];
    }

    #endregion
}