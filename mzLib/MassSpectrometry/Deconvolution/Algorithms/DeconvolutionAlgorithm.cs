﻿using System.Collections.Generic;
using MzLibUtil;

namespace MassSpectrometry
{
    /// <summary>
    /// Parent class defining minimum requirement to be used <see cref="Deconvoluter"/> 
    /// </summary>
    public abstract class DeconvolutionAlgorithm
    {
        protected readonly AverageResidue AverageResidueModel;
        protected readonly DeconvolutionParameters DeconvolutionParameters;

        /// <summary>
        /// Constructor for deconvolution algorithms, nothing should be added to child constructors
        /// </summary>
        /// <param name="deconParameters">parameters to use for deconvolution</param>
        /// <exception cref="MzLibException">thrown in deconvolution parameters did not instantiate fields required by algorithm</exception>
        protected DeconvolutionAlgorithm(DeconvolutionParameters deconParameters)
        {
            DeconvolutionParameters = deconParameters;
            AverageResidueModel = deconParameters.AverageResidueModel;
        }

        /// <summary>
        /// Deconvolutes a mass spectrum
        /// </summary>
        /// <param name="spectrum">spectrum to be deconvoluted</param>
        /// <param name="range">Range of peaks to deconvolute</param>
        /// <returns></returns>
        internal abstract IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range);
    }
}
