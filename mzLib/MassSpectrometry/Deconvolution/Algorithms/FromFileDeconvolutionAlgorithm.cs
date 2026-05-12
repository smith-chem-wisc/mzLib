using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using MzLibUtil;

namespace MassSpectrometry
{
    [ExcludeFromCodeCoverage]
    internal class FromFileDeconvolutionDeconvolutionAlgorithm : DeconvolutionAlgorithm
    {
        internal FromFileDeconvolutionDeconvolutionAlgorithm(DeconvolutionParameters deconParameters) : base(deconParameters)
        {

        }

        internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range = null)
        {
            var deconParams = DeconvolutionParameters as FromFileDeconvolutionParameters ?? throw new MzLibException("Deconvolution params and algorithm do not match");

            var rangeWithRt = range as MzRtRange ?? throw new MzLibException("Deconvolution algorithm only accepts MzRtRange as input range");

            // use range mz and rt to filter results stored in teh parameters. 

            // return those that match the range.

            throw new NotImplementedException();
        }

    }
}
