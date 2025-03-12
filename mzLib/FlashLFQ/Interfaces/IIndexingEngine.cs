using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ.Interfaces
{
    /// <summary>
    /// IIndexingEngine defines the behaviour needed to efficiently read in and index peaks from a mass spectrometric data file
    /// in such a way that they can quickly and efficiently be accessed by calling GetIndexedPeak
    /// </summary>
    public interface IIndexingEngine
    {
        public Ms1ScanInfo[] Ms1ScanInfoArray { get; }
        public SpectraFileInfo SpectraFile { get; }
        public bool IndexPeaks(bool silent);
        public void ClearIndex();
        public void SerializeIndex();
        public void DeserializeIndex();
        public IIndexedPeak GetIndexedPeak(double theorMass, int zeroBasedScanIndex, Tolerance tolerance, int chargeState);
    }
}
