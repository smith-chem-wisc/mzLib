using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class IsoDecDeconvolutionParameters : DeconvolutionParameters
    {
        public override DeconvolutionType DeconvolutionType { get; protected set; } =
            DeconvolutionType.IsoDecDeconvolution;

        public IsoDecDeconvolutionParameters(Polarity polarity = Polarity.Positive) 
            : base(1, 100, polarity)
        {

        }
    }
}
