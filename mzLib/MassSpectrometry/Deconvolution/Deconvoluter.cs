using System.Collections.Generic;
using System.Linq;
using Chemistry;
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
            foreach (var isotopicEnvelope in Deconvolute(scan.MassSpectrum, deconvolutionParameters, rangeToGetPeaksFrom)) 
                yield return isotopicEnvelope;
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
            DeconvolutionAlgorithm deconAlgorithm = deconvolutionParameters.DeconvolutionType switch
            {
                DeconvolutionType.ClassicDeconvolution => new ClassicDeconvolutionAlgorithm(deconvolutionParameters),
                DeconvolutionType.ExampleNewDeconvolutionTemplate => new ExampleNewDeconvolutionAlgorithmTemplate(deconvolutionParameters),
                DeconvolutionType.IsoDecDeconvolution => new IsoDecAlgorithm(deconvolutionParameters),
                _ => throw new MzLibException("DeconvolutionType not yet supported")
            };

            // Short circuit deconvolution if it is called on a neutral mass spectrum
            if (spectrum is NeutralMassSpectrum newt)
            {
                for (int i = 0; i < newt.XArray.Length; i++)
                {
                    // skip this peak if it's outside the range of interest (e.g. if we're only interested in deconvoluting a small m/z range)
                    if (!rangeToGetPeaksFrom.Contains(newt.XArray[i].ToMz(newt.Charges[i])))
                        continue; 
                    yield return new IsotopicEnvelope(newt.XArray[i], newt.YArray[i], newt.Charges[i]);
                }
            }
            else
            {
                foreach (var isotopicEnvelope in deconAlgorithm.Deconvolute(spectrum, rangeToGetPeaksFrom).ToList()) 
                    yield return isotopicEnvelope;
            }

            return deconAlgorithm.Deconvolute(spectrum, rangeToGetPeaksFrom);
        }
    }
}
