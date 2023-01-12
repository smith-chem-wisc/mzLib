using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;

namespace MzLibSpectralAveraging
{
    public static class MassSpectrometryExtensions
    {
        /// <summary>
        /// Average an enumerable of MzSpectrum objects
        /// </summary>
        /// <param name="spectraToAverage">Spectra to average</param>
        /// <param name="parameters">Options for how to average spectra</param>
        /// <returns>Averaged MzSpectrum object</returns>
        public static MzSpectrum AverageSpectra(this IEnumerable<MzSpectrum> spectraToAverage,
            SpectralAveragingParameters parameters)
        {
            var xArrays = spectraToAverage.Select(p => p.XArray).ToArray();
            var yArrays = spectraToAverage.Select(p => p.YArray.SubArray(0, p.YArray.Length)).ToArray();

            var xyJagged = SpectralAveraging.AverageSpectra(xArrays, yArrays, parameters);

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
    }
}
