using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Threading.Tasks;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using SpectralAveraging;

namespace Test.AveragingTests
{
    /// <summary>
    /// Verbatim copies of the averaging code as it was in mzLib 1.0.594, before it was rewritten to stop allocating per
    /// bin. The equivalence tests run the current code and these side by side and require bit-identical output. Only the
    /// class names differ; SpectraNormalization and SpectralWeighting were not changed and are called directly.
    /// </summary>
    [ExcludeFromCodeCoverage]
    internal static class LegacySpectraAveraging
    {
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
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    var peaksFromBin = bins[binIncidences[i]];

                    peaksFromBin = LegacyOutlierRejection.RejectOutliers(peaksFromBin, parameters);
                    if (!peaksFromBin.Any()) continue;

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
            {
                var spectraInBin = bin.Value.Select(p => p.SpectraId).ToArray();
                var representativeBin = bin.Value.First();
                for (var i = 0; i < numSpectra; i++)
                    if (!spectraInBin.Contains(i))
                        bin.Value.Add(new BinnedPeak(representativeBin.Bin, bin.Value.Average(p => p.Mz), 0, i));
            }

            return bins;
        }

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
    }

    /// <summary>
    /// Verbatim copy of OutlierRejection as it was in mzLib 1.0.594; see <see cref="LegacySpectraAveraging"/>.
    /// </summary>
    [ExcludeFromCodeCoverage]
    internal static class LegacyOutlierRejection
    {
        public static double[] RejectOutliers(double[] valuesToCheckForRejection, SpectralAveragingParameters parameters)
        {
            var trimmedMzValues = valuesToCheckForRejection;
            switch (parameters.OutlierRejectionType)
            {
                case OutlierRejectionType.NoRejection:
                    break;

                case OutlierRejectionType.MinMaxClipping:
                    trimmedMzValues = MinMaxClipping(valuesToCheckForRejection);
                    break;

                case OutlierRejectionType.PercentileClipping:
                    trimmedMzValues = PercentileClipping(valuesToCheckForRejection, parameters.Percentile);
                    break;

                case OutlierRejectionType.SigmaClipping:
                    trimmedMzValues = SigmaClipping(valuesToCheckForRejection, parameters.MinSigmaValue,
                        parameters.MaxSigmaValue);
                    break;

                case OutlierRejectionType.WinsorizedSigmaClipping:
                    trimmedMzValues = WinsorizedSigmaClipping(valuesToCheckForRejection, parameters.MinSigmaValue,
                        parameters.MaxSigmaValue);
                    break;

                case OutlierRejectionType.AveragedSigmaClipping:
                    trimmedMzValues = AveragedSigmaClipping(valuesToCheckForRejection, parameters.MinSigmaValue,
                        parameters.MaxSigmaValue);
                    break;

                case OutlierRejectionType.BelowThresholdRejection:
                    trimmedMzValues = BelowThresholdRejection(valuesToCheckForRejection);
                    break;
            }

            return trimmedMzValues;
        }

        internal static List<BinnedPeak> RejectOutliers(List<BinnedPeak> peaks, SpectralAveragingParameters parameters)
        {
            var unrejected = RejectOutliers(peaks.Select(p => p.Intensity).ToArray(), parameters);
            for (var i = 0; i < peaks.Count; i++)
                if (!unrejected.Contains(peaks[i].Intensity))
                {
                    peaks.RemoveAt(i);
                    i--;
                }

            return peaks;
        }

        private static double[] MinMaxClipping(double[] initialValues)
        {
            var max = initialValues.Max();
            var min = initialValues.Min();

            return initialValues.Where(p => p < max && p > min).ToArray();
        }

        private static double[] PercentileClipping(double[] initialValues, double percentile)
        {
            var trim = (1 - percentile) / 2;
            var highPercentile = 1 - trim;
            var median = initialValues.Median();
            var highCutoff = median * (1 + highPercentile);
            var lowCutoff = median * (1 - highPercentile);
            // round to 4-6 decimal places
            return initialValues.Where(p => highCutoff > p && p > lowCutoff).ToArray();
        }

        private static double[] SigmaClipping(double[] initialValues, double sValueMin, double sValueMax)
        {
            var values = initialValues.ToList();
            int n;
            do
            {
                var median = values.Median();
                var standardDeviation = values.StandardDeviation();
                n = 0;
                for (var i = 0; i < values.Count; i++)
                    if (SigmaClipping(values[i], median, standardDeviation, sValueMin, sValueMax))
                    {
                        values.RemoveAt(i);
                        n++;
                        i--;
                    }
            } while (n > 0);

            var val = values.ToArray();
            return val;
        }

        private static double[] WinsorizedSigmaClipping(double[] initialValues, double sValueMin, double sValueMax)
        {
            var values = initialValues.ToList();
            int n;
            const double iterationLimitForHuberLoop = 0.00005;
            do
            {
                if (!values.Any())
                    break;
                var median = values.Median();
                var standardDeviation = values.Where(p => p != 0).StandardDeviation();
                double[] toProcess = new double[values.Count];
                values.CopyTo(toProcess, 0);
                double winsorizedStandardDeviation;

                do // calculates a new median and standard deviation based on the values to do sigma clipping with (Huber loop)
                {
                    var medianLeftBound = median - 1.5 * standardDeviation; // magic number 1.5: optimized parameter per the documentation this method was adapted from, see link at top of class
                    var medianRightBound = median + 1.5 * standardDeviation;
                    toProcess.Winsorize(medianLeftBound, medianRightBound);
                    median = toProcess.Median();
                    winsorizedStandardDeviation = standardDeviation;
                    standardDeviation = toProcess.StandardDeviation() * 1.134; // magic number 1.134: optimized parameter per the documentation this method was adapted from, see link at top of class
                } while (Math.Abs(standardDeviation - winsorizedStandardDeviation) / winsorizedStandardDeviation >
                         iterationLimitForHuberLoop);

                n = 0;
                for (var i = 0; i < values.Count; i++)
                    if (SigmaClipping(values[i], median, standardDeviation, sValueMin, sValueMax))
                    {
                        values.RemoveAt(i);
                        n++;
                        i--;
                    }
            } while (n > 0 && values.Count > 1); // break loop if nothing was rejected, or only one value remains

            return values.ToArray();
        }

        private static double[] AveragedSigmaClipping(double[] initialValues, double sValueMin, double sValueMax)
        {
            var values = initialValues.ToList();
            var median = initialValues.Median();

            double sum = 0;
            for (int i = 0; i < values.Count; i++)
            {
                sum += Math.Pow((values[i] - median), 2) / median;
            }

            double s = Math.Sqrt(sum / (double)(values.Count - 1));

            int n;
            do
            {
                median = values.Where(p => p != 0).Median();
                double sigma = s * Math.Sqrt(median);

                n = 0;
                for (var i = 0; i < values.Count; i++)
                {
                    if (SigmaClipping(values[i], median, sigma, sValueMin, sValueMax))
                    {
                        values.RemoveAt(i);
                        n++;
                        i--;
                    }
                }
            } while (n > 0);

            return values.ToArray();
        }

        private static double[] BelowThresholdRejection(double[] initialValues, double cutoffValue = 0.2)
        {
            var scanCount = initialValues.Length;
            if (initialValues.Count(p => p != 0) <= scanCount * cutoffValue) initialValues = new double[scanCount];
            return initialValues;
        }

        private static void Winsorize(this double[] initialValues, double medianLeftBound, double medianRightBound)
        {
            for (var i = 0; i < initialValues.Length; i++)
            {
                if (initialValues[i] < medianLeftBound)
                {
                    initialValues[i] = medianLeftBound;
                }

                if (initialValues[i] > medianRightBound)
                {
                    initialValues[i] = medianRightBound;
                }
            }
        }

        private static bool SigmaClipping(double value, double median, double standardDeviation, double sValueMin,
            double sValueMax)
        {
            if ((median - value) / standardDeviation > sValueMin)
                return true;
            if ((value - median) / standardDeviation > sValueMax)
                return true;
            return false;
        }
    }

    /// <summary>
    /// Verbatim copy of SpectraFileAveraging.AverageDdaScans and its helper as they were in mzLib 1.0.594; see
    /// <see cref="LegacySpectraAveraging"/>. The spectra themselves are averaged by the current code.
    /// </summary>
    [ExcludeFromCodeCoverage]
    internal static class LegacySpectraFileAveraging
    {
        internal static MsDataScan[] AverageDdaScans(List<MsDataScan> scans, SpectralAveragingParameters parameters)
        {
            List<MsDataScan> averagedScans = new();
            List<MsDataScan> scansToProcess = new();

            int representativeScanMs1Index = parameters.NumberOfScansToAverage % 2 == 0 ? // central scan
                parameters.NumberOfScansToAverage / 2 - 1 : parameters.NumberOfScansToAverage / 2;

            // iterate through all MS1 scans and average them
            foreach (var cvGroup in scans
                .Where(p => p.MsnOrder == 1)
                .OrderBy(p => p.OneBasedScanNumber)
                .GroupBy(p => p.CompensationVoltage))
            {
                foreach (var scan in cvGroup)
                {
                    scansToProcess.Add(scan);
                    // average with new scan from iteration, then remove first scan from list
                    if (scansToProcess.Count != parameters.NumberOfScansToAverage) continue;

                    MsDataScan centralScan = scansToProcess[representativeScanMs1Index];
                    var averagedSpectrum = scansToProcess.AverageSpectra(parameters);
                    var averagedScan = GetAveragedDataScanFromAveragedSpectrum(averagedSpectrum, centralScan);

                    averagedScans.Add(averagedScan);
                    scansToProcess.RemoveAt(0);
                }
                scansToProcess.Clear();
            }

            // add all scans that did not get averaged
            // this includes the MS1 scans from start and end of file and all MS2+ scans
            foreach (var unaveragedScan in scans.Where(original =>
                         !averagedScans.Select(avg => avg.OneBasedScanNumber).Contains(original.OneBasedScanNumber)))
                averagedScans.Add(unaveragedScan);

            return averagedScans.OrderBy(p => p.OneBasedScanNumber).ToArray();
        }

        private static MsDataScan GetAveragedDataScanFromAveragedSpectrum(MzSpectrum averagedSpectrum,
            MsDataScan centralScan)
        {
            MsDataScan averagedScan = new(averagedSpectrum,
                centralScan.OneBasedScanNumber,
                1,
                centralScan.IsCentroid,
                centralScan.Polarity,
                centralScan.RetentionTime,
                averagedSpectrum.Range,
                centralScan.ScanFilter,
                centralScan.MzAnalyzer,
                averagedSpectrum.SumOfAllY,
                centralScan.InjectionTime,
                centralScan.NoiseData,
                centralScan.NativeId,
                centralScan.SelectedIonMZ,
                centralScan.SelectedIonChargeStateGuess,
                centralScan.SelectedIonIntensity,
                centralScan.IsolationMz,
                centralScan.IsolationWidth,
                centralScan.DissociationType,
                centralScan.OneBasedPrecursorScanNumber,
                centralScan.SelectedIonMonoisotopicGuessMz,
                centralScan.HcdEnergy,
                centralScan.ScanDescription,
                centralScan.CompensationVoltage
                );


            var newNativeId =
                averagedScan.NativeId.Replace(averagedScan.NativeId.Split("=").Last(), centralScan.OneBasedScanNumber.ToString());
            averagedScan.SetNativeID(newNativeId);
            return averagedScan;
        }
    }
}
