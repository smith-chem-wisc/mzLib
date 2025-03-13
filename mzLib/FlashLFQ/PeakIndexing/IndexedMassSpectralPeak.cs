using System;
using System.Collections.Generic;

namespace FlashLFQ
{
    [Serializable]
    public class IndexedMassSpectralPeak : ISingleScanDatum, IIndexedPeak
    {
        public int ZeroBasedMs1ScanIndex { get; init; }
        public double Mz { get; init; }
        public double RetentionTime { get; init; }
        public double Intensity { get; protected set; }
        // ISingleScanDatum properties
        public double RelativeSeparationValue => RetentionTime;
        public int ZeroBasedScanIndex => ZeroBasedMs1ScanIndex;

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
            return HashCode.Combine(Mz, ZeroBasedMs1ScanIndex);
        }

        public override string ToString()
        {
            return Mz.ToString("F3") + "; " + ZeroBasedMs1ScanIndex;
        }
    }

    public class IndexedIonMobilityPeak : IndexedMassSpectralPeak, IIndexedPeak
    {
        public HashSet<int> IonMobilityValues { get; init; }
        public int ApexIonMobilityValue { get; init; }

        public IndexedIonMobilityPeak(double mz, double intensity, int zeroBasedMs1ScanIndex, double retentionTime, HashSet<int> ionMobilityValues, int apexIonMobilityValue) 
            : base(mz, intensity, zeroBasedMs1ScanIndex, retentionTime)
           
        {
            this.Mz = mz;
            this.ZeroBasedMs1ScanIndex = zeroBasedMs1ScanIndex;
            this.RetentionTime = retentionTime;
            this.Intensity = intensity;
            this.IonMobilityValues = ionMobilityValues;
            this.ApexIonMobilityValue = apexIonMobilityValue;
        }
    }
}