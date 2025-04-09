using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ.Interfaces
{
    /// <summary>
    /// FlashLFQ requires indexing engines to IndexPeaks, serialize and deserialize indexed peaks,
    /// and retrieve specific peaks (m/z, index pair) or xics (m/z at and around a given retention time
    /// </summary>
    public interface IFlashLfqIndexingEngine
    {
        public ScanInfo[] ScanInfoArray { get; }
        public IIndexedMzPeak GetIndexedPeak(double mz, int zeroBasedScanIndex, PpmTolerance tolerance);
        public List<IIndexedMzPeak> GetXic(double mz, double retentionTime, PpmTolerance ppmTolerance, int missedScansAllowed, double maxPeakHalfWidth = int.MaxValue);
        public void ClearIndex();
        public void SerializeIndex();
        public void DeserializeIndex();
    }
}
