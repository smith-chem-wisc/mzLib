namespace FlashLFQ
{
    public class IsotopicEnvelope
    {
        public readonly IndexedMassSpectralPeak IndexedPeak;
        public readonly int ChargeState;
        public readonly double RetentionTime;

        public IsotopicEnvelope(IndexedMassSpectralPeak monoisotopicPeak, int chargeState, double intensity, double retentionTime)
        {
            IndexedPeak = monoisotopicPeak;
            ChargeState = chargeState;
            Intensity = intensity / chargeState;
            RetentionTime = retentionTime;
        }

        public double Intensity { get; private set; }

        public void Normalize(double normalizationFactor)
        {
            Intensity *= normalizationFactor;
        }
    }
}