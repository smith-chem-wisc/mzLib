using System;

namespace FlashLFQ
{
    [Serializable]
    public class IndexedMassSpectralPeak
    {
        public readonly int ZeroBasedMs1ScanIndex;
        public readonly double Mz;
        public readonly double RetentionTime;
        public readonly double Intensity;

        public IndexedMassSpectralPeak(double mz, double intensity, int zeroBasedMs1ScanIndex, double retentionTime)
        {
            this.Mz = mz;
            this.ZeroBasedMs1ScanIndex = zeroBasedMs1ScanIndex;
            this.RetentionTime = retentionTime;
            this.Intensity = intensity;
        }

        /// <summary>
        /// Creates a deep copy of an IndexedMassSpectralPeak
        /// </summary>
        /// <param name="peak"> IndexedMassSpectralPeak to be copied </param>
        public IndexedMassSpectralPeak(IndexedMassSpectralPeak peak)
        {
            this.Mz = peak.Mz;
            this.ZeroBasedMs1ScanIndex = peak.ZeroBasedMs1ScanIndex;
            this.RetentionTime = peak.RetentionTime;
            this.Intensity = peak.Intensity;
        }

        public override bool Equals(object obj)
        {
            var otherPeak = (IndexedMassSpectralPeak)obj;

            return otherPeak != null
                && otherPeak.Mz == this.Mz
                && otherPeak.ZeroBasedMs1ScanIndex == this.ZeroBasedMs1ScanIndex;
        }

        public override int GetHashCode()
        {
            return Mz.GetHashCode();
        }

        public override string ToString()
        {
            return Mz.ToString("F3") + "; " + ZeroBasedMs1ScanIndex;
        }
    }
}