using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;

namespace MassSpectrometry
{
    public class ThrashDeconvolutionAlgorithm : DeconvolutionAlgorithm
    {
        public ThrashDeconvolutionAlgorithm(DeconvolutionParameters deconParameters) : base(deconParameters)
        {

        }

        public override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            var envelopes = new List<IsotopicEnvelope>();




            return envelopes;
        }
    }
}
