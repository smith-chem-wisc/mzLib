namespace FlashLFQ
{
    public class IsotopicEnvelope
    {
        public readonly IndexedMassSpectralPeak indexedPeak;
        public readonly int chargeState;
        public readonly double retentionTime;
        public double intensity { get; private set; }

        public IsotopicEnvelope(IndexedMassSpectralPeak monoisotopicPeak, int chargeState, double intensity, double retentionTime)
        {
            this.indexedPeak = monoisotopicPeak;
            this.chargeState = chargeState;
            this.intensity = intensity / chargeState;
            this.retentionTime = retentionTime;
        }

        public void Normalize(double normalizationFactor)
        {
            intensity *= normalizationFactor;
        }
    }
}