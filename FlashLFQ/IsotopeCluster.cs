namespace FlashLFQ
{
    public class IsotopeCluster
    {
        #region Public Fields

        public readonly IndexedMassSpectralPeak peakWithScan;
        public readonly int chargeState;
        public double isotopeClusterIntensity;

        #endregion Public Fields

        #region Public Constructors

        public IsotopeCluster(IndexedMassSpectralPeak monoisotopicPeak, int chargeState, double intensity)
        {
            this.peakWithScan = monoisotopicPeak;
            this.chargeState = chargeState;
            this.isotopeClusterIntensity = intensity / chargeState;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            return isotopeClusterIntensity + "; " + peakWithScan.retentionTime;
        }

        #endregion Public Methods
    }
}