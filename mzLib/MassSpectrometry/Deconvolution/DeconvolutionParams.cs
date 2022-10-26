using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;

namespace MassSpectrometry
{
    /// <summary>
    /// Class for hosting all Deconvolution Parameters, regardless of type of deconvolution being performed
    /// </summary>
    public class DeconvolutionParams
    {
        public MzRange? Range { get; set; }
        public int MinAssumedChargeState { get; set; }
        public int MaxAssumedChargeState { get; set; }
        public double DeconvolutionTolerancePpm { get; set; }
        public double? IntensityRatioLimit { get; set; }

        /// <summary>
        /// Constructor should initialize all fields that are used by every deconvolution algorithm
        /// </summary>
        /// <param name="minCharge"></param>
        /// <param name="maxCharge"></param>
        /// <param name="deconPpm"></param>
        public DeconvolutionParams(int minCharge, int maxCharge, double deconPpm)
        {
            MinAssumedChargeState = minCharge;
            MaxAssumedChargeState = maxCharge;
            DeconvolutionTolerancePpm = deconPpm;
        }
    }
}
