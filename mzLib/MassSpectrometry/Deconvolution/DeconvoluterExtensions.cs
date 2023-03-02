using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;

namespace MassSpectrometry
{
    /// <summary>
    /// Class to provide extensions for other objects besides MsDataScan to be deconvoluted
    /// Methods can be generated for each deconvolution type, passing parameters that would otherwise be found in the MsDataScan
    /// </summary>
    public static class DeconvoluterExtensions
    {
        /// <summary>
        /// Deconvolutes a MzSpectrum object
        /// </summary>
        /// <param name="deconvoluter">performs deconvolution</param>
        /// <param name="spectrum">spectrum to deconvolute</param>
        /// <param name="range">mz range of returned peaks</param>
        /// <returns></returns>
        /// <exception cref="MzLibException"></exception>
        public static IEnumerable<IsotopicEnvelope> ClassicDeconvoluteMzSpectra(this Deconvoluter deconvoluter,
            MzSpectrum spectrum, MzRange range = null)
        {
            if (deconvoluter.DeconvolutionType != DeconvolutionTypes.ClassicDeconvolution)
            {
                throw new MzLibException("Deconvoluter is not of correct type for this extension method");
            }
            else
            {
                range ??= spectrum.Range;
                return deconvoluter.DeconvolutionAlgorithm.Deconvolute(spectrum, range);
            }
        }

        public static IEnumerable<IsotopicEnvelope> SpectralDeconvoluteMzSpectra(this Deconvoluter deconvoluter,
            MzSpectrum spectrum, MzRange range)
        {
            if (deconvoluter.DeconvolutionType != DeconvolutionTypes.SpectralDeconvolution)
            {
                throw new MzLibException("Deconvoluter is not of correct type for this extension method");
            }
            else
            {
                ((SpectralDeconvolutionParameters)deconvoluter.DeconvolutionParameters).ScanRange = range;
                return deconvoluter.DeconvolutionAlgorithm.Deconvolute(spectrum);

            }
        }
    }
}
