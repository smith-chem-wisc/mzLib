#nullable enable
using Chemistry;

namespace MassSpectrometry
{
    public abstract class DeconvolutionParameters
    {
        public abstract DeconvolutionType DeconvolutionType { get; protected set; }
        public int MinAssumedChargeState { get; set; }
        public int MaxAssumedChargeState { get; set; }
        public Polarity Polarity { get; set; }
        public AverageResidue AverageResidueModel { get; set; }

        /// <summary>
        /// The expected spacing between isotope peaks in Daltons.
        /// For real (target) deconvolution this is <see cref="Constants.C13MinusC12"/> (~1.003355 Da).
        /// For decoy deconvolution this is a physically impossible value such as 0.9444 Da.
        /// </summary>
        public double ExpectedIsotopeSpacing { get; set; }

        protected DeconvolutionParameters(int minCharge, int maxCharge,
            Polarity polarity = Polarity.Positive,
            AverageResidue? averageResidueModel = null,
            double expectedIsotopeSpacing = Constants.C13MinusC12)
        {
            MinAssumedChargeState = minCharge;
            MaxAssumedChargeState = maxCharge;
            Polarity = polarity;
            AverageResidueModel = averageResidueModel ?? new Averagine();
            ExpectedIsotopeSpacing = expectedIsotopeSpacing;
        }

        /// <summary>
        /// Returns a version of these parameters configured for decoy deconvolution.
        /// The returned object should be lazily constructed and cached — calling this
        /// multiple times should return the same instance.
        /// Subclasses that do not support decoy deconvolution should return <c>null</c>.
        /// </summary>
        /// <returns>
        /// A <see cref="DeconvolutionParameters"/> configured for decoy spacing,
        /// or <c>null</c> if this parameter type does not support decoy deconvolution.
        /// </returns>
        public abstract DeconvolutionParameters? ToDecoyParameters();
    }
}