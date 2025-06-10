#nullable enable

namespace MassSpectrometry
{
    /// <summary>
    /// Classic Deconvolution Required Parameters
    /// </summary>
    public class ClassicDeconvolutionParameters : DeconvolutionParameters
    {
        public override DeconvolutionType DeconvolutionType { get; protected set; } = DeconvolutionType.ClassicDeconvolution;
        public double DeconvolutionTolerancePpm { get; set; }
        public double IntensityRatioLimit { get; set; }

        /// <summary>
        /// Construct Classic deconvolution parameters
        /// </summary>
        public ClassicDeconvolutionParameters(int minCharge, int maxCharge, double deconPpm, double intensityRatio, Polarity polarity = Polarity.Positive, AverageResidue? averageResidueModel = null)
            : base(minCharge, maxCharge, polarity, averageResidueModel)
        {
            IntensityRatioLimit = intensityRatio;
            DeconvolutionTolerancePpm = deconPpm;
        }
    }
}
