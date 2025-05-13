#nullable enable
namespace MassSpectrometry
{
    /// <summary>
    /// Class for hosting deconvolution parameters common to all methods
    /// </summary>
    public abstract class DeconvolutionParameters
    {
        public abstract DeconvolutionType DeconvolutionType { get; protected set; }
        public int MinAssumedChargeState { get; set; }
        public int MaxAssumedChargeState { get; set; }
        public Polarity Polarity { get; set; }
        public AverageResidue AverageResidueModel { get; set; }

        /// <summary>
        /// Constructor should initialize all fields that are used by every deconvolution algorithm
        /// </summary>
        public DeconvolutionParameters(int minCharge, int maxCharge, Polarity polarity = Polarity.Positive, AverageResidue? averageResidueModel = null)
        {
            MinAssumedChargeState = minCharge;
            MaxAssumedChargeState = maxCharge;
            Polarity = polarity;
            AverageResidueModel = averageResidueModel ?? new Averagine(); // Default to Averagine
        }
    }
}
    

