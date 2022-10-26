using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;

namespace MassSpectrometry
{
    public class AlexDeconv : DeconvolutionAlgorithm
    {
        public AlexDeconv(DeconvolutionParams deconParams) : base(deconParams)
        {

        }

        public override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum)
        {
            throw new NotImplementedException();
        }

        protected override bool CheckAlgorithmParameterCompatibility()
        {
            throw new NotImplementedException();
        }
    }
}
