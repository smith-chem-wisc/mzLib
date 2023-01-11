using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using MzLibUtil;

namespace MzLibSpectralAveraging
{
    public static class SpectralAveraging
    {
        /// <summary>
        /// Average a group of spectra in the jagged array format
        /// </summary>
        /// <param name="xArrays"> x values of spectra</param>
        /// <param name="yArrays">y values of spectra</param>
        /// <param name="parameters">Options for how to perform averaging</param>
        /// <returns>Averaged MzSpectrum object</returns>
        /// <exception cref="NotImplementedException">If merging type has not been implemented</exception>
        public static double[][] AverageSpectra(double[][] xArrays, double[][] yArrays, SpectralAveragingParameters parameters)
        {
            switch (parameters.SpectraMergingType)
            {
                case SpectraMergingType.MzBinning:
                     return MzBinning(xArrays, yArrays, parameters);

                default:
                    throw new NotImplementedException("Spectrum Merging Type Not Yet Implemented");
            }
        }

        #region Extension Methods

        /// <summary>
        /// Average an enumerable of MzSpectrum objects
        /// </summary>
        /// <param name="spectraToAverage">Spectra to average</param>
        /// <param name="parameters">Options for how to average spectra</param>
        /// <returns>Averaged MzSpectrum object</returns>
        public static MzSpectrum AverageSpectra(this IEnumerable<MzSpectrum> spectraToAverage,
            SpectralAveragingParameters parameters)
        {
            var toAverage = spectraToAverage as MzSpectrum[] ?? spectraToAverage.ToArray();
            var xArrays = toAverage.Select(p => p.XArray).ToArray();
            var yArrays = toAverage.Select(p => p.YArray.SubArray(0, p.YArray.Length)).ToArray();

            //var xyJagged = AverageSpectra(xArrays, yArrays, options);

            var xyJagged = SpectrumBinning(xArrays, yArrays, parameters);

            return new MzSpectrum(xyJagged[0], xyJagged[1], true);
        }

        /// <summary>
        /// Average an enumerable of MsDataScans objects
        /// </summary>
        /// <param name="scansToAverage">Scans to average</param>
        /// <param name="parameters">Options for how to average scans</param>
        /// <returns>Averaged MzSpectrum object</returns>
        public static MzSpectrum AverageSpectra(this IEnumerable<MsDataScan> scansToAverage, SpectralAveragingParameters parameters)
        {
            return scansToAverage.Select(p => p.MassSpectrum).AverageSpectra(parameters);
        }

        #endregion

        private static double[][] MzBinning(double[][] xArrays, double[][] yArrays, SpectralAveragingParameters parameters)
        {
            PixelStack pixelStack = new(xArrays, yArrays, parameters.BinSize);
            if (parameters.PerformNormalization) pixelStack.PerformNormalization();
            SpectralWeighting.CalculateSpectraWeights(pixelStack, parameters);
            OutlierRejection.RejectOutliers(pixelStack, parameters);
            return pixelStack.MergeSpectra(parameters);
        }



        public static double[][] SpectrumBinning(double[][] xArrays, double[][] yArrays,
            SpectralAveragingParameters parameters)
        {
            // get tics 
            double[] tics = yArrays.Select(p => p.Sum()).ToArray();
            double averageTic = tics.Average();

            // normalize spectra
            if (parameters.NormalizationType != NormalizationType.NoNormalization)
            {
                SpectraNormalization.NormalizeSpectra(yArrays, parameters.NormalizationType);
            }

            // get bins
            var bins = GetBins(xArrays, yArrays, parameters.BinSize);

            // get weights
            var weights = SpectralWeighting.WeightEvenly(yArrays.Length);

            // reject outliers and average bins
            List <(double mz, double intensity)> averagedPeaks = new();
            foreach (var bin in bins)
            {
                bins[bin.Key] = OutlierRejection.RejectOutliers(bin.Value, parameters);
                averagedPeaks.Add(AverageBin(bin.Value, weights));
            }

            // return averaged
            return new[]
            {
                averagedPeaks.OrderBy(p => p.mz).Select(p => p.mz).ToArray(),
                averagedPeaks.Select(p => p.intensity).ToArray()
            };

        }

        private static Dictionary<int, List<BinnedPeak>> GetBins(double[][] xArrays, double[][] yArrays,
            double binSize)
        {
            var numSpectra = xArrays.Length;

            // calculate the number of bins to be utilized
            double min = 100000;
            for (int i = 0; i < numSpectra; i++)
            {
                min = Math.Min(xArrays[i][0], min);
            }

            Dictionary<int, List<BinnedPeak>> bins = new();
            for (int i = 0; i < numSpectra; i++)
            {
                for (int j = 0; j < xArrays[i].Length; j++)
                {
                    if (yArrays[i][j] == 0) continue;

                    int binIndex = (int)Math.Floor((xArrays[i][j] - min) / binSize);
                    var binValue = new BinnedPeak(binIndex, xArrays[i][j], yArrays[i][j], i);
                    if (!bins.ContainsKey(binIndex))
                    {
                        bins.Add(binIndex, new List<BinnedPeak>() {binValue});
                    }
                    else
                    {
                        bins[binIndex].Add(binValue);
                    }
                }
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
            var intensity = numerator / denominator;
            var mz = peaksInBin.Select(p => p.Mz).Average();

            return (mz, intensity);
        }


        /// <summary>
        /// Merges spectra into a two dimensional array of (m/z, int) values based upon their bin 
        /// </summary>
        /// <param name="scans">scans to be combined</param>
        /// <returns>MSDataScan with merged values</returns>
        public static double[][] SpectrumBinning(double[][] xArrays, double[][] yArrays, double[] totalIonCurrents, double binSize, int numSpectra,
            SpectralAveragingParameters parameters)
        {
            // normalize if selected
            if (parameters.PerformNormalization)
            {
                double averageIonCurrent = totalIonCurrents.Average();
                for (int i = 0; i < xArrays.Length; i++)
                {
                    for (int j = 0; j < xArrays[i].Length; j++)
                    {
                        yArrays[i][j] /= totalIonCurrents[i] * averageIonCurrent;
                    }
                }
            }

            // calculate the bins to be utilized
            double min = 100000;
            double max = 0;
            for (int i = 0; i < numSpectra; i++)
            {
                min = Math.Min(xArrays[i][0], min);
                max = Math.Max(xArrays[i].Max(), max);
            }
            int numberOfBins = (int)Math.Ceiling((max - min) * (1 / binSize));

            double[][] xValuesBin = new double[numberOfBins][];
            double[][] yValuesBin = new double[numberOfBins][];
            // go through each scan and place each (m/z, int) from the spectra into a jagged array
            for (int i = 0; i < numSpectra; i++)
            {
                for (int j = 0; j < xArrays[i].Length; j++)
                {
                    int binIndex = (int)Math.Floor((xArrays[i][j] - min) / binSize);
                    if (xValuesBin[binIndex] == null)
                    {
                        xValuesBin[binIndex] = new double[numSpectra];
                        yValuesBin[binIndex] = new double[numSpectra];
                    }
                    xValuesBin[binIndex][i] = xArrays[i][j];
                    yValuesBin[binIndex][i] = yArrays[i][j];
                }
            }

            xValuesBin = xValuesBin.Where(p => p != null).ToArray();
            yValuesBin = yValuesBin.Where(p => p != null).ToArray();

            // average the remaining arrays to create the composite spectrum
            // this will clipping and averaging for y values as indicated in the settings
            double[] xArray = new double[xValuesBin.Length];
            double[] yArray = new double[yValuesBin.Length];
            // target for optimization here
            for (int i = 0; i < yValuesBin.Length; i++)
            {
                // linq is probably slow 
                xArray[i] = xValuesBin[i].Where(p => p != 0).Average();
                yArray[i] = ProcessSingleMzArray(yValuesBin[i].OrderBy(p => p).ToArray(), parameters);
            }

            return new double[][] { xArray, yArray };
        }

        /// <summary>
        /// Main Engine of this Binning method, processes a single array of intesnity values for a single mz and returns their average
        /// </summary>
        /// <param name="intInitialArray"></param>
        /// <returns></returns>
        public static double ProcessSingleMzArray(double[] intInitialArray, SpectralAveragingParameters parameters)
        {
            double average;
            double[] weights;
            double[] trimmedArray;

            if (intInitialArray.Where(p => p != 0).Count() <= 1)
                return 0;

            trimmedArray = OutlierRejection.RejectOutliers(intInitialArray, parameters);
            if (trimmedArray.Where(p => p != 0).Count() <= 1)
                return 0;
            weights = Enumerable.Repeat(1.0, intInitialArray.Length).ToArray();
            average = MergePeakValuesToAverage(trimmedArray, weights);
            return average;
        }

        /// <summary>
        /// Calculates the weighted average value for each m/z point passed to it
        /// </summary>
        /// <param name="intValues">array of mz values to evaluate </param>
        /// <param name="weights">relative weight assigned to each of the mz values</param>
        /// <returns></returns>
        public static double MergePeakValuesToAverage(double[] intValues, double[] weights)
        {
            double numerator = 0;
            for (int i = 0; i < intValues.Count(); i++)
            {
                numerator += intValues[i] * weights[i];
            }
            double average = numerator / weights.Sum();
            return average;
        }
    }
}
