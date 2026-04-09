#nullable enable
using Chemistry;

namespace MassSpectrometry
{
    public class ClassicDeconvolutionParameters : DeconvolutionParameters
    {
        public override DeconvolutionType DeconvolutionType { get; protected set; }
            = DeconvolutionType.ClassicDeconvolution;
        public double DeconvolutionTolerancePpm { get; set; }
        public double IntensityRatioLimit { get; set; }

        public ClassicDeconvolutionParameters(int minCharge, int maxCharge,
            double deconPpm, double intensityRatio,
            Polarity polarity = Polarity.Positive,
            AverageResidue? averageResidueModel = null,
            double expectedIsotopeSpacing = Constants.C13MinusC12)
            : base(minCharge, maxCharge, polarity, averageResidueModel, expectedIsotopeSpacing)
        {
            IntensityRatioLimit = intensityRatio;
            DeconvolutionTolerancePpm = deconPpm;
        }

        private DeconvolutionParameters? _decoyParams = null;
        public override DeconvolutionParameters ToDecoyParameters() =>
            _decoyParams ??= new ClassicDeconvolutionParameters(
                MinAssumedChargeState,
                MaxAssumedChargeState,
                DeconvolutionTolerancePpm,
                IntensityRatioLimit,
                Polarity,
                averageResidueModel: new DecoyAveragine(AverageResidueModel,
                    DecoyAveragine.DefaultDecoyIsotopeSpacing),
                expectedIsotopeSpacing: DecoyAveragine.DefaultDecoyIsotopeSpacing);
    }
}