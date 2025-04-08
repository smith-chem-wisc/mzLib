using MzLibUtil;
using System;
using System.Collections.Generic;
using MassSpectrometry;

namespace FlashLFQ.Interfaces
{
    /// <summary>
    /// IIndexingEngine defines the behaviour needed to efficiently read in and index peaks from a mass spectrometric data file
    /// in such a way that they can quickly and efficiently be accessed by calling GetIndexedPeak
    /// </summary>
    public interface IIndexingEngine
    {
        public ScanInfo[] ScanInfoArray { get; }
        public bool IndexPeaks(MsDataScan[] scanArray);
        public IIndexedMzPeak GetIndexedPeak(double mz, int zeroBasedScanIndex, PpmTolerance tolerance);
        public List<IIndexedMzPeak> GetXic(double mz, double retentionTime, PpmTolerance ppmTolerance, int missedScansAllowed, double maxPeakHalfWidth = int.MaxValue);
        public List<IIndexedMzPeak> GetXic(double mz, int zeroBasedStartIndex, PpmTolerance ppmTolerance, int missedScansAllowed, double maxPeakHalfWidth = int.MaxValue);
    }
}
