#nullable enable
using Chemistry;

namespace MassSpectrometry
{
    /// <summary>
    /// Class for hosting deconvolution parameters common to all methods.
    /// </summary>
    public abstract class DeconvolutionParameters
    {
        public abstract DeconvolutionType DeconvolutionType { get; protected set; }
        public int MinAssumedChargeState { get; set; }
        public int MaxAssumedChargeState { get; set; }
        public Polarity Polarity { get; set; }
        public AverageResidue AverageResidueModel { get; set; }

        /// <summary>
        /// The isotope spacing (in Da) used during deconvolution.
        /// The default value — <see cref="Constants.C13MinusC12"/> ≈ 1.003355 Da —
        /// corresponds to the mass difference between ¹³C and ¹²C and is correct for
        /// all standard deconvolution runs.
        ///
        /// For target-decoy FDR estimation, this is set to 0.9444 Da by
        /// <see cref="DeconvolutionDecoyGenerator.MakeDecoyParameters"/> to generate
        /// envelopes with a physically impossible isotope spacing (the canonical
        /// value from OpenMS FLASHDeconv). All other parameters are unchanged, so
        /// decoy envelopes are subject to the same quality filters as targets.
        ///
        /// Do not change this value manually unless you are implementing a new
        /// decoy strategy. Setting it to anything other than the default will produce
        /// incorrect deconvolution results.
        /// </summary>
        public double DecoyIsotopeDistance { get; set; } = Constants.C13MinusC12;

        /// <summary>
        /// Marks this parameters instance as belonging to a decoy deconvolution pass.
        /// Informational only — no algorithm logic branches on this value.
        /// Set to <c>true</c> by <see cref="DeconvolutionDecoyGenerator.MakeDecoyParameters"/>.
        /// </summary>
        public bool IsDecoyRun { get; set; } = false;

        /// <summary>
        /// Constructor should initialize all fields that are used by every deconvolution algorithm.
        /// </summary>
        protected DeconvolutionParameters(int minCharge, int maxCharge, Polarity polarity = Polarity.Positive, AverageResidue? averageResidueModel = null)
        {
            MinAssumedChargeState = minCharge;
            MaxAssumedChargeState = maxCharge;
            Polarity = polarity;
            AverageResidueModel = averageResidueModel ?? new Averagine(); // Default to Averagine
        }

        /// <summary>
        /// Creates a shallow clone of this instance. Used by
        /// <see cref="DeconvolutionDecoyGenerator.MakeDecoyParameters"/> to produce a
        /// decoy parameters object of the same concrete type without reflection.
        ///
        /// Subclasses that hold mutable reference-type state beyond the base class fields
        /// should override this method to perform a deeper copy of those fields.
        /// The base implementation is sufficient for all current mzLib parameter types
        /// because all non-array fields are value types or immutable.
        /// </summary>
        internal virtual DeconvolutionParameters ShallowClone()
            => (DeconvolutionParameters)MemberwiseClone();
    }
}
