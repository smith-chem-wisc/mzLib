using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ.Interfaces
{
    public interface IIndexingEngine
    {
        public Ms1ScanInfo[] Ms1ScanInfoArray { get; }
        public bool IndexPeaks(SpectraFileInfo fileInfo, bool silent);
        public void ClearIndex();
        public void SerializeIndex(SpectraFileInfo file);
        public void DeserializeIndex(SpectraFileInfo file);
        public IIndexedPeak GetIndexedPeak(double mz, int zeroBasedScanIndex, Tolerance tolerance);
    }
}
