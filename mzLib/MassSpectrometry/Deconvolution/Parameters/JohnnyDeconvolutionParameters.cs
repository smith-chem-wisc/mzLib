using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class JohnnyDeconvolutionParameters : DeconvolutionParameters
    {
        public override DeconvolutionType DeconvolutionType { get; protected set; }

        public JohnnyDeconvolutionParameters(int minCharge, int maxCharge, Polarity polarity = Polarity.Positive) 
            : base(minCharge, maxCharge, polarity)
        {
        }
    }
}
