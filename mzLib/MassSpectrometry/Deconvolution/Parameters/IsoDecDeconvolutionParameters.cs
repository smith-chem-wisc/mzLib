using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class IsoDecDeconvolutionParameters(Polarity polarity = Polarity.Positive)
        : DeconvolutionParameters(1, 100, polarity)
    {
        public override DeconvolutionType DeconvolutionType { get; protected set; } = DeconvolutionType.IsoDecDeconvolution;
    }
}
