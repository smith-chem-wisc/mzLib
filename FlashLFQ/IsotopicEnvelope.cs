namespace FlashLFQ
{
    public class IsotopicEnvelope
    {
        #region Public Fields

        public readonly IndexedMassSpectralPeak indexedPeak;
        public readonly int chargeState;
        public readonly double retentionTime;
        public double intensity { get; private set; }

        #endregion Public Fields

        #region Public Constructors

        public IsotopicEnvelope(IndexedMassSpectralPeak monoisotopicPeak, int chargeState, double intensity, double retentionTime)
        {
            this.indexedPeak = monoisotopicPeak;
            this.chargeState = chargeState;
            this.intensity = intensity / chargeState;
            this.retentionTime = retentionTime;
        }
        
        #endregion Public Constructors

        #region Public Methods

        public void Normalize(double normalizationFactor)
        {
            intensity *= normalizationFactor;
        }

        #endregion Public Methods
    }
}