using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MzLibSpectralAveraging
{
    public static class SpectralAveraging
    {
        /// <summary>
        /// Average a group of spectra in the jagged array format
        /// </summary>
        /// <param name="xArrays"> x values of spectra</param>
        /// <param name="yArrays">y values of spectra</param>
        /// <param name="options">Options for how to perform averaging</param>
        /// <returns>Averaged MzSpectrum object</returns>
        /// <exception cref="NotImplementedException">If merging type has not been implemented</exception>
        public static double[][] AverageSpectra(double[][] xArrays, double[][] yArrays, SpectralAveragingOptions options)
        {
            switch (options.SpectraMergingType)
            {
                case SpectraMergingType.MzBinning:
                     return MzBinning(xArrays, yArrays, options);

                default:
                    throw new NotImplementedException("Spectrum Merging Type Not Yet Implemented");
            }
        }

        #region Extension Methods

        /// <summary>
        /// Average an enumerable of MzSpectrum objects
        /// </summary>
        /// <param name="spectraToAverage">Spectra to average</param>
        /// <param name="options">Options for how to average spectra</param>
        /// <returns>Averaged MzSpectrum object</returns>
        public static MzSpectrum AverageSpectra(this IEnumerable<MzSpectrum> spectraToAverage,
            SpectralAveragingOptions options)
        {
            var xArrays = spectraToAverage.Select(p => p.XArray).ToArray();
            var yArrays = spectraToAverage.Select(p => p.YArray).ToArray();

            var xyJagged = AverageSpectra(xArrays, yArrays, options);
            return new MzSpectrum(xyJagged[0], xyJagged[1], true);
        }

        /// <summary>
        /// Average an enumerable of MsDataScans objects
        /// </summary>
        /// <param name="scansToAverage">Scans to average</param>
        /// <param name="options">Options for how to average scans</param>
        /// <returns>Averaged MzSpectrum object</returns>
        public static MzSpectrum AverageSpectra(this IEnumerable<MsDataScan> scansToAverage, SpectralAveragingOptions options)
        {
            return scansToAverage.Select(p => p.MassSpectrum).AverageSpectra(options);
        }

        #endregion

        private static double[][] MzBinning(double[][] xArrays, double[][] yArrays, SpectralAveragingOptions options)
        {
            BinnedSpectra binnedSpectra = new(xArrays, yArrays, options.BinSize);
            if (options.PerformNormalization) binnedSpectra.PerformNormalization();
            SpectralWeighting.CalculateSpectraWeights(binnedSpectra, options);
            OutlierRejection.RejectOutliers(binnedSpectra, options);
            return binnedSpectra.MergeSpectra(options);
        }
    }
}
