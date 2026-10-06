using FlashLFQ.Interfaces;
using MassSpectrometry;
using MzLibUtil;
using NetSerializer;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;

namespace FlashLFQ
{
    public class PeakIndexingEngine : IndexingEngine<IndexedMassSpectralPeak>, IFlashLfqIndexingEngine
    {
        private readonly Serializer _serializer;

        /// <summary>
        /// The index SerializeIndex set aside in memory instead of writing to disk. Null when the
        /// index was written to disk, or has been restored by DeserializeIndex.
        /// </summary>
        private List<IndexedMassSpectralPeak>[]? _setAsideIndex;

        /// <summary>
        /// SerializeIndex keeps an index in memory while the machine's memory load, including that
        /// index, stays at or below this fraction of the memory available to the process.
        /// </summary>
        internal const double MaxMemoryLoadFractionToKeepIndex = 0.5;

        /// <summary>
        /// The size of one IndexedMassSpectralPeak on a 64-bit runtime: a 16-byte header and four 4-byte fields
        /// </summary>
        private const long BytesPerIndexedPeak = 32;

        public SpectraFileInfo SpectraFile { get; private set; }
        internal PeakIndexingEngine()
        {
            var messageTypes = new List<Type>
            {
                typeof(List<IndexedMassSpectralPeak>[]), typeof(List<IndexedMassSpectralPeak>),
                typeof(IndexedMassSpectralPeak)
            };
            _serializer = new Serializer(messageTypes);
        }

        /// <summary>
        /// This factory method returns an IndexingEngine instance where the peaks in all MS1 scans have been indexed. 
        /// This method ignores MS2 scans when indexing
        /// </summary>
        public static PeakIndexingEngine? InitializeIndexingEngine(SpectraFileInfo file)
        {
            // read spectra file
            string fileName = file.FullFilePathWithExtension;
            var reader = MsDataFileReader.GetDataFile(fileName);
            reader.LoadAllStaticData();

            var peakIndexingEngine = InitializeIndexingEngine(reader);
            if(peakIndexingEngine != null) peakIndexingEngine.SpectraFile = file;
            return peakIndexingEngine;
        }

        /// <summary>
        /// This factory method returns an IndexingEngine instance where the peaks in all MS1 scans have been indexed. 
        /// This method ignores MS2 scans when indexing
        /// </summary>
        public static PeakIndexingEngine? InitializeIndexingEngine(MsDataFile dataFile)
        {
            var scanArray = dataFile.GetMS1Scans()
                .Where(i => i != null && i.MsnOrder == 1)
                .OrderBy(i => i.OneBasedScanNumber)
                .ToArray();
            return InitializeIndexingEngine(scanArray);
        }

        /// <summary>
        /// Read in all spectral peaks from the scanArray, index the peaks based on mass and retention time, 
        /// and store them in a jagged array of Lists containing all peaks within a particular mass range
        /// </summary>
        /// <param name="scanArray">An array of raw data scans</param>
        public static PeakIndexingEngine? InitializeIndexingEngine(MsDataScan[] scanArray)
        {
            PeakIndexingEngine newEngine = new();
            if (newEngine.IndexPeaks(scanArray))
                return newEngine;
            return null;
        }

        public void ClearIndex()
        {
            bool releasesMemory = IndexedPeaks != null && !ReferenceEquals(IndexedPeaks, _setAsideIndex);
            IndexedPeaks = null;
            // Collecting is pointless when the index is still held in memory to be restored later
            if (releasesMemory)
                GC.Collect();
        }

        /// <summary>
        /// Sets the index aside until DeserializeIndex restores it. The index stays in memory when
        /// there is room for it (see <see cref="MaxMemoryLoadFractionToKeepIndex"/>); otherwise it is
        /// written to a .ind file in the same directory as the spectra file.
        /// </summary>
        public void SerializeIndex()
        {
            GCMemoryInfo memoryInfo = GC.GetGCMemoryInfo();
            if (IndexFitsInMemory(EstimateIndexBytes(), memoryInfo.MemoryLoadBytes, memoryInfo.TotalAvailableMemoryBytes))
            {
                _setAsideIndex = IndexedPeaks;
                return;
            }
            WriteIndexToDisk();
        }

        /// <summary>
        /// Writes the index to a .ind file in the same directory as the spectra file
        /// </summary>
        internal void WriteIndexToDisk()
        {
            string dir = Path.GetDirectoryName(SpectraFile.FullFilePathWithExtension);
            string indexPath = Path.Combine(dir, SpectraFile.FilenameWithoutExtension + ".ind");

            using (var indexFile = File.Create(indexPath))
            {
                _serializer.Serialize(indexFile, IndexedPeaks);
            }
        }

        /// <summary>
        /// Restores the index set aside by SerializeIndex, from memory or from its .ind file.
        /// </summary>
        public void DeserializeIndex()
        {
            if (_setAsideIndex != null)
            {
                IndexedPeaks = _setAsideIndex;
                _setAsideIndex = null;
                return;
            }

            string dir = Path.GetDirectoryName(SpectraFile.FullFilePathWithExtension);
            string indexPath = Path.Combine(dir, SpectraFile.FilenameWithoutExtension + ".ind");

            using (var indexFile = File.OpenRead(indexPath))
            {
                IndexedPeaks = (List<IndexedMassSpectralPeak>[])_serializer.Deserialize(indexFile);
            }

            File.Delete(indexPath);
        }

        /// <summary>
        /// True when keeping an index of <paramref name="indexBytes"/> in memory leaves the memory
        /// load at or below <see cref="MaxMemoryLoadFractionToKeepIndex"/> of the available memory.
        /// An unknown available memory (zero, before the first garbage collection) never fits.
        /// </summary>
        internal static bool IndexFitsInMemory(long indexBytes, long memoryLoadBytes, long totalAvailableMemoryBytes)
        {
            return totalAvailableMemoryBytes > 0
                && memoryLoadBytes + indexBytes <= totalAvailableMemoryBytes * MaxMemoryLoadFractionToKeepIndex;
        }

        /// <summary>
        /// The approximate managed size of the current index, in bytes
        /// </summary>
        internal long EstimateIndexBytes()
        {
            if (IndexedPeaks == null)
                return 0;
            long peaks = 0;
            long references = IndexedPeaks.Length;
            foreach (var bin in IndexedPeaks)
            {
                if (bin == null) continue;
                peaks += bin.Count;
                references += bin.Capacity;
            }
            return peaks * BytesPerIndexedPeak + references * IntPtr.Size;
        }

        /// <summary>  
        /// Prune the index engine to remove any unnecessary data or entries for the better memory usage.  
        /// </summary>  
        public void PruneIndex(List<float> targetMz)
        {
            PpmTolerance ppmTolerance = new PpmTolerance(10); // Default tolerance, can be adjusted as needed
            if (IndexedPeaks == null || targetMz == null || !targetMz.Any())
                return;
            List<IndexedMassSpectralPeak>[] indexedPeaks = new List<IndexedMassSpectralPeak>[IndexedPeaks.Length];
            var maxIndex = IndexedPeaks.Length - 1;

            foreach (var mz in targetMz)
            {
                int ceilingMz = (int)Math.Ceiling(ppmTolerance.GetMaximumValue(mz) * BinsPerDalton);
                int floorMz = (int)Math.Floor(ppmTolerance.GetMinimumValue(mz) * BinsPerDalton);
                if (ceilingMz > maxIndex || floorMz > maxIndex)
                {
                    // If the mz is out of bounds, skip it
                    continue;
                }

                for (int i = floorMz; i <= ceilingMz; i++)
                {
                    if (indexedPeaks[i] == null)
                    {
                        indexedPeaks[i] = IndexedPeaks[i];
                    }

                }
            }
            IndexedPeaks = indexedPeaks;
        }
    }
}