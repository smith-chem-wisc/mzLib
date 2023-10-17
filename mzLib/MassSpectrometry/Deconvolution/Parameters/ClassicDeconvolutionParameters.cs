using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;

namespace MassSpectrometry
{
    /// <summary>
    /// Classic Deconvolution Required Parameters
    /// </summary>
    public class ClassicDeconvolutionParameters : DeconvolutionParameters
    {
        public double DeconvolutionTolerancePpm { get; set; }
        public double IntensityRatioLimit { get; set; }

        /// <summary>
        /// Construct Classic deconvolution parameters
        /// </summary>
        /// <param name="minCharge"></param>
        /// <param name="maxCharge"></param>
        /// <param name="deconPpm"></param>
        /// <param name="intensityRatio"></param>
        /// <param name="range">Isolation range of the scan to be deconvoluted</param>
        public ClassicDeconvolutionParameters(int minCharge, int maxCharge, double deconPpm, double intensityRatio, Polarity polarity = Polarity.Positive)
            : base(minCharge, maxCharge, polarity)
        {
            IntensityRatioLimit = intensityRatio;
            DeconvolutionTolerancePpm = deconPpm;
        }
    }
}
