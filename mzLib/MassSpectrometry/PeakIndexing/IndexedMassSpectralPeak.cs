using System;

namespace MassSpectrometry
{
    [Serializable]
    public class IndexedMassSpectralPeak : IIndexedPeak
    {
        public int ZeroBasedScanIndex { get; init; }
        public float Mz { get; init; }
        public float M => Mz;
        public float RetentionTime { get; init; }
        public float Intensity { get; init; }
        public IndexedMassSpectralPeak(float mz, float intensity, int zeroBasedMs1ScanIndex, float retentionTime)
        {
            this.Mz = mz;
            this.ZeroBasedScanIndex = zeroBasedMs1ScanIndex;
            this.RetentionTime = retentionTime;
            this.Intensity = intensity;
        }
        public override bool Equals(object obj)
        {
            return Equals((IndexedMassSpectralPeak)obj);
        }

        public bool Equals(IIndexedPeak other)
        {
            return Equals((IndexedMassSpectralPeak)other);
        }

        public bool Equals(IndexedMassSpectralPeak other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return other.ZeroBasedScanIndex == this.ZeroBasedScanIndex
                   && Math.Abs(other.Mz - this.Mz) < 1e-9;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(Mz, ZeroBasedScanIndex);
        }

        public override string ToString()
        {
            return Mz.ToString("F3") + "; " + ZeroBasedScanIndex;
        }
    }
}