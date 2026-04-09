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

        // Thread-safe lazy caching of decoy parameters using double-checked locking.
        // The ??= operator alone is not atomic: two threads can both observe null,
        // each construct a new instance, and the last writer wins — meaning the first
        // thread's caller gets an object that isn't the one ultimately cached.
        // The volatile field ensures cross-thread visibility of writes, the outer null
        // check provides a fast path (no lock contention once cached), and the inner
        // ??= inside the lock guarantees only one thread ever constructs the instance.
        private volatile DeconvolutionParameters? _decoyParams = null;
        private readonly object _decoyParamsLock = new();
        public override DeconvolutionParameters ToDecoyParameters()
        {
            if (_decoyParams is not null)
                return _decoyParams;

            lock (_decoyParamsLock)
            {
                return _decoyParams ??= new ClassicDeconvolutionParameters(
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
    }
}