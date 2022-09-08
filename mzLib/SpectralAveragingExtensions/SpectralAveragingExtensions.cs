using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using SpectralAveraging;

namespace SpectralAveragingExtensions
{
    public static class SpectralAveragingExtensions
    {
        #region Averaging

        /// <summary>
        /// Average an enumerable of MzSpectrum objects
        /// </summary>
        /// <param name="spectraToAverage">Spectra to average</param>
        /// <param name="options">Options for how to average spectra</param>
        /// <returns></returns>
        public static MzSpectrum CombineSpectra(this List<MzSpectrum> spectraToAverage,
            SpectralAveragingOptions options)
        {
            double[][] compositeSpectraValues = SpectralMerging.CombineSpectra(spectraToAverage.Select(p => p.XArray).ToArray(),
                spectraToAverage.Select(p => p.YArray).ToArray(), spectraToAverage.Select(p => p.SumOfAllY).ToArray(),
                spectraToAverage.Count(), options);
            return new MzSpectrum(compositeSpectraValues[0], compositeSpectraValues[1], true);
        }

        /// <summary>
        /// Average an enumerable of MsDataScans objects
        /// </summary>
        /// <param name="scansToAverage">Scans to average</param>
        /// <param name="options">Options for how to average scans</param>
        /// <returns></returns>
        public static MzSpectrum CombineSpectra(this List<MsDataScan> scansToAverage, SpectralAveragingOptions options)
        {
            return scansToAverage.Select(p => p.MassSpectrum).ToList().CombineSpectra(options);
        }

        #endregion

        #region Normalization

        /// <summary>
        /// Normalize many spectra in MzSpectrum format
        /// </summary>
        /// <param name="spectraToNormalize">Spectra to normalize</param>
        /// <param name="multiplyByAverageTic">true will make it a relative normalization, false will make it an absolute normalization</param>
        public static void NormalizeSpectrumToTic(this List<MzSpectrum> spectraToNormalize, bool multiplyByAverageTic)
        {
            double averageTic = spectraToNormalize.Select(p => p.SumOfAllY).Average();
            for (int i = 0; i < spectraToNormalize.Count(); i++)
            {
                if (multiplyByAverageTic)
                {
                    SpectrumNormalization.NormalizeSpectrumToTic(spectraToNormalize[i].YArray,
                        spectraToNormalize[i].SumOfAllY, averageTic);
                }
                else
                {
                    SpectrumNormalization.NormalizeSpectrumToTic(spectraToNormalize[i].YArray,
                        spectraToNormalize[i].SumOfAllY);
                }
            }
        }

        /// <summary>
        /// Normalize many spectra in MsDataScan format
        /// </summary>
        /// <param name="scansToNormalize">Scans to normalize</param>
        /// <param name="multiplyByAverageTic">true will make it a relative normalization, false will make it an absolute normalization</param>
        public static void NormalizeSpectrumToTic(this List<MsDataScan> scansToNormalize, bool multiplyByAverageTic)
        {
            scansToNormalize.Select(p => p.MassSpectrum).ToList().NormalizeSpectrumToTic(multiplyByAverageTic);
        }

        /// <summary>
        /// Normalize single MzSpectrum
        /// Total of all intensities will be near 1
        /// </summary>
        /// <param name="spectrumToNormalize"></param>
        public static void NormalizeSpectrumToTic(this MzSpectrum spectrumToNormalize)
        {
            SpectrumNormalization.NormalizeSpectrumToTic(spectrumToNormalize.YArray, spectrumToNormalize.SumOfAllY);
        }

        /// <summary>
        /// Normalize single MsDataScan
        /// Sum of all intensities will near one
        /// </summary>
        /// <param name="scanToNormalize"></param>
        public static void NormalizeSpectrumToTic(this MsDataScan scanToNormalize)
        {
            scanToNormalize.MassSpectrum.NormalizeSpectrumToTic();
        }

        #endregion

    }
}
