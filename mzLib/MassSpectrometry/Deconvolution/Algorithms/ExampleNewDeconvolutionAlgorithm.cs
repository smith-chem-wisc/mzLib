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
    public class ExampleNewDeconvolutionAlgorithm : DeconvolutionAlgorithm
    {
        public ExampleNewDeconvolutionAlgorithm(DeconvolutionParameters deconParameters) : base(deconParameters)
        {

        }

        public override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum)
        {
            var deconParams = DeconvolutionParameters as ExampleNewDeconvolutionParameters ?? throw new MzLibException("Deconvolution params and algorithm do not match");
            throw new NotImplementedException();
        }

    }
}
