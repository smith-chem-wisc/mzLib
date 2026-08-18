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
        public SpectraFileInfo SpectraFile { get; }
        /// <summary>
        /// A generic method for finding the closest peak with a specified m/z and in a specified scan. Returns null if no peaks within tolerance are found.
        /// </summary>
        /// <param name="mz"> the m/z of the peak to be searched for </param>
        /// <param name="zeroBasedScanIndex"> the zero based index of the scan where the peak is to be found </param>
        /// <param name="charge"> an optional parameter used only for IIndexedMass and massIndexingEngine; must be null for mz peak indexing </param>
        public IIndexedPeak GetIndexedPeak(double mz, int zeroBasedScanIndex, Tolerance tolerance, int? charge = null);
        /// <summary>
        /// A generic method of peak tracing across the retention time. Finds peaks with a given mz that occur on either side of a given
        /// retention time. Peak searching iterates backwards through the scans until the peak 
        /// is no longer observed (i.e., is absent in more scans than allowed, as defined by the
        /// missedScansAllowed parameter. Missed scans don't have to be sequential. The same procedure
        /// is then repeated in the forward direction.
        /// </summary>
        /// <param name="mz"> the m/z of the peak to be searched for </param>
        /// <param name="retentionTime"> the retention time where peak searching will begin </param>
        /// <param name="missedScansAllowed"> the number of successive missed scans allowed before the xic is terminated </param>
        /// <param name="maxPeakHalfWidth"> the maximum distance from the apex RT of the XIC to both start RT and end RT </param>
        /// <param name="matchedPeaks"> the dictionary that stores all the peaks already matched to an xic </param>
        /// <param name="charge"> an optional parameter used only for IIndexedMass and massIndexingEngine; must be null for mz peak indexing </param>
        /// <returns> A list of IIndexedPeak objects, ordered by retention time </returns>
        public List<IIndexedPeak> GetXic(double mz, double retentionTime, Tolerance ppmTolerance, int missedScansAllowed, double maxPeakHalfWidth = int.MaxValue, int? charge = null, Dictionary<IIndexedPeak, ExtractedIonChromatogram> matchedPeaks = null);
        /// <summary>
        /// Clear the indexed peaks and the jagged array of indexed peaks to free up memory
        /// </summary>
        public void ClearIndex();
        /// <summary>
        /// Write the indexed peaks to a file. File will be written to the same directory as the data file, as determined by the SpectraFileInfo property
        /// </summary>
        public void SerializeIndex();
        /// <summary>
        /// Reads the indexed peaks from a file. File will be read from the same directory as the data file, as determined by the SpectraFileInfo property
        /// </summary>
        public void DeserializeIndex();

        /// <summary>
        /// Prune the index engine to remove any unnecessary data or entries for the better memory usage.
        /// </summary>
        public void PruneIndex(List<float> targetMass);
    }
}
