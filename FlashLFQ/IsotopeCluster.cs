namespace FlashLFQ
{
    public class IsotopeCluster
    {
        #region Public Fields

        public readonly IndexedMassSpectralPeak indexedPeak;
        public readonly int chargeState;
        public readonly double retentionTime;
        public double isotopeClusterIntensity;

        #endregion Public Fields

        #region Public Constructors

        public IsotopeCluster(IndexedMassSpectralPeak monoisotopicPeak, int chargeState, double intensity, double retentionTime)
        {
            this.indexedPeak = monoisotopicPeak;
            this.chargeState = chargeState;
            this.isotopeClusterIntensity = intensity / chargeState;
            this.retentionTime = retentionTime;
        }

        #endregion Public Constructors
    }
}