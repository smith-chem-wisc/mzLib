using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Easy.Common.Extensions;
using MzLibUtil;

namespace MassSpectrometry
{
    public enum DeconvolutionType
    {
        ClassicDeconvolution,
        ExampleNewDeconvolution,
    }

    /// <summary>
    /// Context class for all deconvolution
    /// </summary>
    public class Deconvoluter
    {
        public DeconvolutionAlgorithm DeconvolutionAlgorithm { get; private set; }
        public DeconvolutionType DeconvolutionType { get; }
        public DeconvolutionParameters DeconvolutionParameters { get; }

        public Deconvoluter(DeconvolutionType deconType, DeconvolutionParameters deconParameters)
        {
            DeconvolutionParameters = deconParameters;
            DeconvolutionType = deconType;
            ConstructDeconvolutionAlgorithm(deconParameters);
        }

        /// <summary>
        /// Deconvolute a MsDataScan
        /// </summary>
        /// <param name="scan">scan to deconvolute</param>
        /// <param name="rangeToGetPeaksFrom">Range of peaks to deconvolute, if null, will deconvolute entire spectra</param>
        /// <returns></returns>
        public IEnumerable<IsotopicEnvelope> Deconvolute(MsDataScan scan, MzRange rangeToGetPeaksFrom = null)
        {
            rangeToGetPeaksFrom ??= scan.MassSpectrum.Range;

            // set deconvolution parameters that are only present in the MsDataScan
            switch (DeconvolutionType)
            {
                case DeconvolutionType.ClassicDeconvolution:
                    break;

                case DeconvolutionType.ExampleNewDeconvolution:
                    break;
            }

            return DeconvolutionAlgorithm.Deconvolute(scan.MassSpectrum, rangeToGetPeaksFrom);
        }

 

        /// <summary>
        /// Constructs the relevant deconvolution algorithm
        /// </summary>
        /// <param name="deconParameters"></param>
        /// <exception cref="MzLibException">if a type of enum is used that is not supported</exception>
        private void ConstructDeconvolutionAlgorithm(DeconvolutionParameters deconParameters)
        {
            // construct algorithm object
            switch (DeconvolutionType)
            {
                case DeconvolutionType.ClassicDeconvolution:
                    DeconvolutionAlgorithm = new ClassicDeconvolutionAlgorithm(deconParameters);
                    break;

                case DeconvolutionType.ExampleNewDeconvolution:
                    DeconvolutionAlgorithm = new ExampleNewDeconvolutionAlgorithm(deconParameters);
                    break;

                default: throw new MzLibException("DeconvolutionType not yet supported");
            }
        }
    }
}
