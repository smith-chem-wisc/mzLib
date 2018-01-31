namespace FlashLFQ
{
    public class IsotopeCluster
    {
        public readonly IndexedMassSpectralPeak peakWithScan;
        public readonly int chargeState;
        public double isotopeClusterIntensity;

        public IsotopeCluster(IndexedMassSpectralPeak monoisotopicPeak, int chargeState, double intensity)
        {
            this.peakWithScan = monoisotopicPeak;
            this.chargeState = chargeState;
            this.isotopeClusterIntensity = intensity / chargeState;
        }

        public override string ToString()
        {
            return isotopeClusterIntensity + "; " + peakWithScan.retentionTime;
        }
    }
}
