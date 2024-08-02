using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class JohnnyDeconvolutionParameters : DeconvolutionParameters
    {
        public override DeconvolutionType DeconvolutionType { get; protected set; } =
            DeconvolutionType.JohnnyDeconvolution;

        public JohnnyDeconvolutionParameters(Polarity polarity = Polarity.Positive) 
            : base(1, 100, polarity)
        {

        }
    }
}
