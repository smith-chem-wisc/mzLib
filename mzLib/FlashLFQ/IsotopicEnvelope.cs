namespace FlashLFQ
{
    /// <summary>
    /// Contains the summed intensities of all isotope peaks detected in a single MS1 scan for a given species.
    /// </summary>
    public class IsotopicEnvelope
    {
        /// <summary>
        /// The most abundant isotopic peak used for peak finding.
        /// </summary>
        public readonly IndexedMassSpectralPeak IndexedPeak;
        public readonly int ChargeState;

        public IsotopicEnvelope(IndexedMassSpectralPeak monoisotopicPeak, int chargeState, double intensity)
        {
            IndexedPeak = monoisotopicPeak;
            ChargeState = chargeState;
            Intensity = intensity / chargeState;
        }

        /// <summary>
        /// The summed intensity of all isotope peaks detected in one MS1 scan. This sum may contain 
        /// imputed intensity values for expected isotopes that weren't observed, but only if the observed 
        /// isotopic distribution was otherwise similar to the expected isotopic distribution.
        /// </summary>
        public double Intensity { get; private set; }

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