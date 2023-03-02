﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Easy.Common.Extensions;
using MassSpectrometry.Deconvolution.Algorithms;
using MzLibUtil;

namespace MassSpectrometry
{
    public enum DeconvolutionTypes
    {
        ClassicDeconvolution,
        SpectralDeconvolution,
    }

    /// <summary>
    /// Context class for all deconvolution
    /// </summary>
    public class Deconvoluter
    {
        public DeconvolutionAlgorithm DeconvolutionAlgorithm { get; private set; }
        public DeconvolutionTypes DeconvolutionType { get; }
        public DeconvolutionParameters DeconvolutionParameters { get; }

        public Deconvoluter(DeconvolutionTypes deconType, DeconvolutionParameters deconParameters)
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
                case DeconvolutionTypes.ClassicDeconvolution:
                    ((ClassicDeconvolutionParameters)DeconvolutionParameters).Range = 
                        new MzRange(scan.IsolationRange.Minimum - 8.5, scan.IsolationRange.Maximum + 8.5);
                    break;

                case DeconvolutionTypes.SpectralDeconvolution:
                    ((SpectralDeconvolutionParameters)DeconvolutionParameters).ScanRange =
                        new MzRange(scan.IsolationRange.Minimum - 8.5, scan.IsolationRange.Maximum + 8.5);
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
                case DeconvolutionTypes.ClassicDeconvolution:
                    DeconvolutionAlgorithm = new ClassicDeconvolutionAlgorithm(deconParameters);
                    break;

                case DeconvolutionTypes.SpectralDeconvolution:
                    DeconvolutionAlgorithm = new SpectralDeconvolutionAlgorithm(deconParameters);
                    break;

                default: throw new MzLibException("DeconvolutionType not yet supported");
            }
        }
    }
}
