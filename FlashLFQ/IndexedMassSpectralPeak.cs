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