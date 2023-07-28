using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    [ExcludeFromCodeCoverage]
    public class ExampleNewDeconvolutionParameters : DeconvolutionParameters
    {
        public ExampleNewDeconvolutionParameters(int minCharge, int maxCharge, Polarity polarity = Polarity.Positive)
            : base(minCharge, maxCharge, polarity)
        {

        }
    }
}
