using MzLibUtil;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    public class MultipleDeconvolutionAlgorithm : DeconvolutionAlgorithm
    {
        public MultipleDeconvolutionAlgorithm(DeconvolutionParameters deconParameters) : base(deconParameters)
        {

        }

        protected internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range = null)
        {
            var deconParams = DeconvolutionParameters as MultipleDeconParameters ?? throw new MzLibException("Deconvolution params and algorithm do not match");
            range ??= spectrum.Range;

            IEnumerable<IsotopicEnvelope> envelopes = [];
            foreach (var param in deconParams.Parameters)
            {
                var result = Deconvoluter.Deconvolute(spectrum, param, range);
                envelopes = envelopes.Concat(result);
            }

            return envelopes;
        }
    }
}
