using System.Collections.Generic;
using MassSpectrometry.Deconvolution;
using MzLibUtil;

namespace MassSpectrometry
{
    public abstract class DeconvolutionAlgorithm
    {
        public readonly AverageResidue AverageResidueModel;
        protected readonly DeconvolutionParameters DeconvolutionParameters;

        /// <summary>
        /// Constructor for deconvolution algorithms, nothing should be added to child constructors
        /// </summary>
        /// <param name="deconParameters">parameters to use for deconvolution</param>
        /// <exception cref="MzLibException">thrown in deconvolution parameters did not instantiate fields required by algorithm</exception>
        protected DeconvolutionAlgorithm(DeconvolutionParameters deconParameters)
        {
            DeconvolutionParameters = deconParameters;
            AverageResidueModel = deconParameters.Polarity == Polarity.Positive ? new Averagine() : new Averatide();
        }

        /// <summary>
        /// Deconvolutes a mass spectrum
        /// </summary>
        /// <param name="spectrum">spectrum to be deconvoluted</param>
        /// <param name="range">Range of peaks to deconvolute</param>
        /// <returns></returns>
        public abstract IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range);
    }
}
