using FlashLFQ.Interfaces;
using System;
using System.Collections.Generic;

namespace FlashLFQ
{
    [Serializable]
    public class IndexedMassSpectralPeak : IIndexedMzPeak
    {
        public int ZeroBasedScanIndex { get; init; }
        public double Mz { get; init; }
        public double RetentionTime { get; init; }
        public double Intensity { get; init; }

        public IndexedMassSpectralPeak(double mz, double intensity, int zeroBasedMs1ScanIndex, double retentionTime)
        {
            this.Mz = mz;
            this.ZeroBasedScanIndex = zeroBasedMs1ScanIndex;
            this.RetentionTime = retentionTime;
            this.Intensity = intensity;
        }

        public override bool Equals(object obj)
        {
            var otherPeak = (IndexedMassSpectralPeak)obj;

            return otherPeak != null
                && otherPeak.Mz == this.Mz
                && otherPeak.ZeroBasedScanIndex == this.ZeroBasedScanIndex;
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