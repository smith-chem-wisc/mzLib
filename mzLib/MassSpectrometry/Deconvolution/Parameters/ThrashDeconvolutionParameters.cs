using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class ThrashDeconvolutionParameters : DeconvolutionParameters
    {
        public double MinMsFeatureToBackgroundRation { get; set; }

        public ThrashDeconvolutionParameters(double minMsFeatureToBackgroundRation = 3)
        {
            MinMsFeatureToBackgroundRation = minMsFeatureToBackgroundRation;
        }
    }
}
