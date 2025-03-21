using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using MzLibUtil;

namespace MassSpectrometry
{
    [ExcludeFromCodeCoverage]
    internal class ExampleNewDeconvolutionAlgorithmTemplate : DeconvolutionAlgorithm
    {
        internal ExampleNewDeconvolutionAlgorithmTemplate(DeconvolutionParameters deconParameters) : base(deconParameters)
        {

        }

        internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range = null)
        {
            var deconParams = DeconvolutionParameters as ExampleNewDeconvolutionParametersTemplate ?? throw new MzLibException("Deconvolution params and algorithm do not match");
            range ??= spectrum.Range;

            throw new NotImplementedException();
        }

    }
}
