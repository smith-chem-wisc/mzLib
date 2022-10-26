using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Easy.Common.Extensions;
using MzLibUtil;

namespace MassSpectrometry
{
    public enum DeconvolutionTypes
    {
        ClassicDeconvolution,
        AlexDeconvolution,
    }

    /// <summary>
    /// Context class for all deconvolution
    /// </summary>
    public class Deconvoluter
    {
        private readonly DeconvolutionAlgorithm deconvolutionAlgorithm;
        private readonly DeconvolutionParams deconvolutionParams;

        public Deconvoluter(DeconvolutionTypes deconType, DeconvolutionParams deconParams)
        {
            // verify parameters are compatible with algorithm
            deconvolutionParams = deconParams;

            // construct algorithm object
            switch (deconType)
            {
                case DeconvolutionTypes.ClassicDeconvolution:
                    deconvolutionAlgorithm = new ClassicDeconv(deconParams);
                    break;

                case DeconvolutionTypes.AlexDeconvolution:
                    deconvolutionAlgorithm = new AlexDeconv(deconParams);
                    break;

                default: throw new MzLibException("DeconvoltionType not yet supported");
            }
        }

        /// <summary>
        /// Deconvolute a MsDataScan
        /// </summary>
        /// <param name="scan"></param>
        /// <returns></returns>
        public IEnumerable<IsotopicEnvelope> Deconvolute(MsDataScan scan)
        {
            var range =
                new MzRange(scan.IsolationRange.Minimum - 8.5, scan.IsolationRange.Maximum + 8.5);
            return Deconvolute(scan.MassSpectrum, range);
        }

        /// <summary>
        /// Deconvolute an MzSpectrum with its isolation range
        /// NOTE: Range may only be required for ClassicDeconv, not sure yet
        /// if so, this should be removed
        /// </summary>
        /// <param name="spectrum"></param>
        /// <param name="range"></param>
        /// <returns></returns>
        public IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            deconvolutionParams.Range = range;
            return deconvolutionAlgorithm.Deconvolute(spectrum);
        }
    }
}
