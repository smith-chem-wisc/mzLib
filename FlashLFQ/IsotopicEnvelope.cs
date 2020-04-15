namespace FlashLFQ
{
    public class IsotopicEnvelope
    {
        public readonly IndexedMassSpectralPeak IndexedPeak;
        public readonly int ChargeState;

        public IsotopicEnvelope(IndexedMassSpectralPeak monoisotopicPeak, int chargeState, double intensity, double corr)
        {
            IndexedPeak = monoisotopicPeak;
            ChargeState = chargeState;
            Intensity = intensity / chargeState;
            PearsonCorrelationToTheoretical = corr;
        }

        public double Intensity { get; private set; }
        public readonly double PearsonCorrelationToTheoretical;

        public void Normalize(double normalizationFactor)
        {
            Intensity *= normalizationFactor;
        }

        public override string ToString()
        {
            return "+" + ChargeState + "|" + Intensity.ToString("F0") + "|" + IndexedPeak.RetentionTime.ToString("F3") + "|" + IndexedPeak.ZeroBasedMs1ScanIndex;
        }

        public override bool Equals(object obj)
        {
            var otherEnv = (IsotopicEnvelope)obj;

            return otherEnv != null
                && otherEnv.ChargeState == this.ChargeState
                && otherEnv.IndexedPeak.Equals(this.IndexedPeak);
        }

        public override int GetHashCode()
        {
            return ChargeState.GetHashCode() + IndexedPeak.GetHashCode();
        }
    }
}