using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;

namespace MassSpectrometry
{
    [ExcludeFromCodeCoverage]
    public class ExampleNewDeconvolutionAlgorithmTemplate : DeconvolutionAlgorithm
    {
        public ExampleNewDeconvolutionAlgorithmTemplate(DeconvolutionParameters deconParameters) : base(deconParameters)
        {

        }

        public override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range = null)
        {
            var deconParams = DeconvolutionParameters as ExampleNewDeconvolutionParametersTemplate ?? throw new MzLibException("Deconvolution params and algorithm do not match");
            range ??= spectrum.Range;

            throw new NotImplementedException();
        }

    }
}
