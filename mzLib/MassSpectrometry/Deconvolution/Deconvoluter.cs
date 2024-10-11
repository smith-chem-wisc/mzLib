using System.Collections.Generic;
using MzLibUtil;

namespace MassSpectrometry
{
    public enum DeconvolutionType
    {
        ClassicDeconvolution,
        ExampleNewDeconvolutionTemplate,
    }

    /// <summary>
    /// Context class for all deconvolution
    /// </summary>
    public static class Deconvoluter
    {
        /// <summary>
        /// Static deconvolution of an MsDataScan that does not require Deconvoluter construction
        /// </summary>
        /// <param name="scan">scan to deconvolute</param>
        /// <param name="deconvolutionParameters">decon parameters to use, also determines type of deconvolution used</param>
        /// <param name="rangeToGetPeaksFrom">Range of peaks to deconvolute, if null, will deconvolute entire spectra</param>
        /// <returns></returns>
        public static IEnumerable<IsotopicEnvelope> Deconvolute(MsDataScan scan,
            DeconvolutionParameters deconvolutionParameters, MzRange rangeToGetPeaksFrom = null)
        {
            // set any specific deconvolution parameters found only in the MsDataScan
            
            return Deconvolute(scan.MassSpectrum, deconvolutionParameters, rangeToGetPeaksFrom);
        }
                    


        /// <summary>
        /// Static deconvolution of an MzSpectrum that does not require Deconvoluter construction
        /// </summary>
        /// <param name="spectrum">spectrum to deconvolute</param>
        /// <param name="deconvolutionParameters">decon parameters to use, also determines type of deconvolution used</param>
        /// <param name="rangeToGetPeaksFrom">Range of peaks to deconvolute, if null, will deconvolute entire spectra</param>
        /// <returns></returns>
        public static IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum,
            DeconvolutionParameters deconvolutionParameters, MzRange rangeToGetPeaksFrom = null)
        {
            rangeToGetPeaksFrom ??= spectrum.Range;

            // set deconvolution algorithm 
            DeconvolutionAlgorithm deconAlgorithm;
            switch (deconvolutionParameters.DeconvolutionType)
            {
                case DeconvolutionType.ClassicDeconvolution:
                    deconAlgorithm = new ClassicDeconvolutionAlgorithm(deconvolutionParameters);
                    break;

                case DeconvolutionType.ExampleNewDeconvolutionTemplate:
                    deconAlgorithm = new ExampleNewDeconvolutionAlgorithmTemplate(deconvolutionParameters);
                    break;

                default: throw new MzLibException("DeconvolutionType not yet supported");
            }

            return deconAlgorithm.Deconvolute(spectrum, rangeToGetPeaksFrom);
        }
    }
}
