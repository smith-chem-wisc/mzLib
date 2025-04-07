using Chemistry;
using Readers;
using MassSpectrometry;
using MzLibUtil;
using NetSerializer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using FlashLFQ.Interfaces;

namespace FlashLFQ
{
    public class PeakIndexingEngine : IIndexingEngine
    {
        private List<IIndexedMzPeak>[] _indexedPeaks;
        private readonly Serializer _serializer;
        private const int BinsPerDalton = 100;
        public Ms1ScanInfo[] ScanInfoArray { get; private set; }
        public SpectraFileInfo SpectraFile { get; init; }

        public PeakIndexingEngine(SpectraFileInfo file)
        {
            var messageTypes = new List<Type>
            {
                typeof(List<IIndexedMzPeak>[]), typeof(List<IIndexedMzPeak>),
                typeof(IIndexedMzPeak)
            };
            _serializer = new Serializer(messageTypes);
            SpectraFile = file;
        }

        public PeakIndexingEngine(MsDataScan[] scans)
        {
            PeakIndexing(scans);
        }

        public bool IndexPeaks(bool silent)
        {
            if (!silent)
            {
                Console.WriteLine("Reading spectra file");
            }

            // read spectra file
            string fileName = SpectraFile.FullFilePathWithExtension;
            var reader = MsDataFileReader.GetDataFile(fileName); 

            reader.LoadAllStaticData();
            // retrieve only the ms1s. 
            MsDataScan[] msDataScans = reader.GetMS1Scans()
                .Where(i => i != null && i.MsnOrder == 1)
                .OrderBy(i => i.OneBasedScanNumber)
                .ToArray(); 
            

            if (msDataScans.All(p => p == null))
            {
                _indexedPeaks = null;
                return false;
            }

            PeakIndexing(msDataScans);


            if (_indexedPeaks == null || _indexedPeaks.Length == 0)
            {
                if (!silent)
                {
                    Console.WriteLine("FlashLFQ Error: The file " + SpectraFile.FilenameWithoutExtension + " contained no MS1 peaks!");
                }

                return false;
            }

            return true;
        }

        /// <summary>
        /// Read in all spectral peaks from scans, index the peaks and store them in a list ordered by m/z
        /// </summary>
        /// <param name="msDataScans">An array of raw data scans</param>
        /// <param name="scanInfo">Outputs a list of scan information for each scan which is needed for FlashLfq
        public void PeakIndexing(MsDataScan[] msDataScans)
        {
            _indexedPeaks = new List<IIndexedMzPeak>[(int)Math.Ceiling(msDataScans.Where(p => p != null
                && p.MassSpectrum.LastX != null).Max(p => p.MassSpectrum.LastX.Value) * BinsPerDalton) + 1];
            ScanInfoArray = new Ms1ScanInfo[msDataScans.Length];

            for (int scanIndex = 0; scanIndex < msDataScans.Length; scanIndex++)
            {
                ScanInfoArray[scanIndex] = new Ms1ScanInfo(msDataScans[scanIndex].OneBasedScanNumber, scanIndex, msDataScans[scanIndex].RetentionTime);

                for (int j = 0; j < msDataScans[scanIndex].MassSpectrum.XArray.Length; j++)
                {
                    int roundedMz = (int)Math.Round(msDataScans[scanIndex].MassSpectrum.XArray[j] * BinsPerDalton, 0);
                    _indexedPeaks[roundedMz] ??= new List<IIndexedMzPeak>();
                    _indexedPeaks[roundedMz].Add(
                        new IndexedMassSpectralPeak(
                            msDataScans[scanIndex].MassSpectrum.XArray[j],
                            msDataScans[scanIndex].MassSpectrum.YArray[j], 
                            scanIndex, 
                            msDataScans[scanIndex].RetentionTime));
                }

            }
        }

        public void ClearIndex()
        {
            _indexedPeaks = null;
            GC.Collect();
        }

        public void SerializeIndex()
        {
            string dir = Path.GetDirectoryName(SpectraFile.FullFilePathWithExtension);
            string indexPath = Path.Combine(dir, SpectraFile.FilenameWithoutExtension + ".ind");

            using (var indexFile = File.Create(indexPath))
            {
                _serializer.Serialize(indexFile, _indexedPeaks);
            }
            ClearIndex();
        }

        public void DeserializeIndex()
        {
            string dir = Path.GetDirectoryName(SpectraFile.FullFilePathWithExtension);
            string indexPath = Path.Combine(dir, SpectraFile.FilenameWithoutExtension + ".ind");

            using (var indexFile = File.OpenRead(indexPath))
            {
                _indexedPeaks = (List<IIndexedMzPeak>[])_serializer.Deserialize(indexFile);
            }

            File.Delete(indexPath);
        }

        public IIndexedMzPeak GetIndexedPeak(double theorMass, int zeroBasedScanIndex, PpmTolerance ppmTolerance, int chargeState) =>
            GetIndexedPeak(theorMass.ToMz(chargeState), zeroBasedScanIndex, ppmTolerance);

        /// <summary>
        /// A generic method for finding the closest peak with a specified m/z and in a specified scan. Returns null if no peaks within tolerance are found.
        /// </summary>
        /// <param name="mz"> the m/z of the peak to be searched for </param>
        /// <param name="zeroBasedScanIndex"> the zero based index of the scan where the peak is to be found </param>
        public IIndexedMzPeak GetIndexedPeak(double mz, int zeroBasedScanIndex, PpmTolerance ppmTolerance)
        {
            var bins = GetBinsInRange(mz, ppmTolerance);
            if (bins.Count == 0) return null;
            List<int> peakIndicesInBins = bins.Select(b => BinarySearchForIndexedPeak(b, zeroBasedScanIndex)).ToList();
            return GetBestPeakFromBins(bins, mz, zeroBasedScanIndex, peakIndicesInBins, ppmTolerance);
        }

        /// <summary>
        /// A generic method of peak tracing across the retention time. Finds peaks with a given mz that occur on either side of a given
        /// retention time. Peak searching iterates backwards through the scans until the peak 
        /// is no longer observed (i.e., is absent in more scans than allowed, as defined by the
        /// missedScansAllowed parameter. Missed scans don't have to be sequential. The same procedure
        /// is then repeated in the forward direction.
        /// </summary>
        /// <param name="mz"> the m/z of the peak to be searched for </param>
        /// <param name="zeroBasedStartIndex"> the scan where peak searching behaviour begins </param>
        /// <param name="missedScansAllowed"> the number of successive missed scans allowed before the xic is terminated </param>
        /// <param name="maxPeakHalfWidth"> the maximum distance from the apex RT of the XIC to both start RT and end RT </param>
        /// <returns></returns>
        public List<IIndexedMzPeak> GetXic(double mz, int zeroBasedStartIndex, PpmTolerance ppmTolerance, int missedScansAllowed, double maxPeakHalfWidth = double.MaxValue)
        {
            List<IIndexedMzPeak> xic = new List<IIndexedMzPeak>();
            var allBins = GetBinsInRange(mz, ppmTolerance);
            if (allBins.Count == 0)
                return xic;

            // For each bin, find + store a pointer to the current index
            int[] peakPointerArray = allBins.Select(b => BinarySearchForIndexedPeak(b, zeroBasedStartIndex)).ToArray();
            IIndexedMzPeak initialPeak = GetBestPeakFromBins(allBins, mz, zeroBasedStartIndex, peakPointerArray, ppmTolerance);

            if (initialPeak != null)
                xic.Add(initialPeak);

            foreach (int direction in new List<int> { -1, 1 })
            {
                //int missedPeaks = initialPeak == null ? 1 : 0;
                int missedPeaks = 0; // In the legacy code, the initial peak was not counted as a missed peak if it wasn't found. This should change in the future, but for now, we will keep it as is.
                int currentZeroBasedScanIndex = zeroBasedStartIndex;
                var pointerArrayCopy = new int[peakPointerArray.Length];
                Array.Copy(peakPointerArray, pointerArrayCopy, peakPointerArray.Length);

                while (missedPeaks < missedScansAllowed)
                {
                    //increment all pointers
                    for (int i = 0; i < pointerArrayCopy.Length; i++)
                    {
                        pointerArrayCopy[i] += direction;
                    }
                    currentZeroBasedScanIndex += direction;

                    // Search for the next peak
                    IIndexedMzPeak nextPeak = GetBestPeakFromBins(allBins, mz, currentZeroBasedScanIndex, pointerArrayCopy, ppmTolerance);

                    // Add the peak to the XIC or increment the missed peaks
                    if (nextPeak == null)
                        missedPeaks++;
                    else
                        xic.Add(nextPeak);
                }
            }

            // Sort the XIC in place
            xic.Sort((x, y) => x.RetentionTime.CompareTo(y.RetentionTime));
            return xic;
        }

        #region Peak finding helper functions

        internal List<List<IIndexedMzPeak>> GetBinsInRange(double mz, PpmTolerance ppmTolerance)
        {
            int ceilingMz = (int)Math.Ceiling(ppmTolerance.GetMaximumValue(mz) * BinsPerDalton);
            int floorMz = (int)Math.Floor(ppmTolerance.GetMinimumValue(mz) * BinsPerDalton);
            List<List<IIndexedMzPeak>> allBins = new();
            for (int j = floorMz; j <= ceilingMz; j++)
            {
                if (j >= _indexedPeaks.Length || _indexedPeaks[j] == null)
                    continue;
                allBins.Add(_indexedPeaks[j]);
            }
            return allBins;
        }

        internal static IIndexedMzPeak GetBestPeakFromBins(List<List<IIndexedMzPeak>> allBins, double mz, int zeroBasedScanIndex, IList<int> peakIndicesInBins, PpmTolerance ppmTolerance)
        {
            IIndexedMzPeak bestPeak = null;
            for (int i = 0; i < allBins.Count; i++)
            {
                var tempPeak = GetPeakFromBin(allBins[i], mz, zeroBasedScanIndex, peakIndicesInBins[i], ppmTolerance);
                if (ppmTolerance.Within(tempPeak.Mz, mz)
                    && tempPeak.ZeroBasedScanIndex == peakIndicesInBins[i]
                    && (bestPeak == null || Math.Abs(tempPeak.Mz - mz) < Math.Abs(bestPeak.Mz - mz)))
                {
                    bestPeak = tempPeak;
                }
            }
            return bestPeak;
        }

        internal static IIndexedMzPeak GetPeakFromBin(List<IIndexedMzPeak> bin, double mz, int zeroBasedScanIndex, int peakIndexInBin, PpmTolerance ppmTolerance)
        {
            IIndexedMzPeak bestPeak = null;
            for (int i = peakIndexInBin; i < bin.Count; i++)
            {
                IIndexedMzPeak peak = bin[i];

                if (peak.ZeroBasedScanIndex > zeroBasedScanIndex)
                {
                    break;
                }

                if (ppmTolerance.Within(peak.Mz, mz)
                    && peak.ZeroBasedScanIndex == zeroBasedScanIndex
                    && (bestPeak == null || Math.Abs(peak.Mz - mz) < Math.Abs(bestPeak.Mz - mz)))
                {
                    bestPeak = peak;
                }
            }
            return bestPeak;
        }

        internal static int BinarySearchForIndexedPeak(List<IIndexedMzPeak> indexedPeaks, int zeroBasedScanIndex)
        {
            int m = 0;
            int l = 0;
            int r = indexedPeaks.Count - 1;

            while (l <= r)
            {
                m = l + ((r - l) / 2);

                if (r - l < 2)
                {
                    break;
                }
                if (indexedPeaks[m].ZeroBasedScanIndex < zeroBasedScanIndex)
                {
                    l = m + 1;
                }
                else
                {
                    r = m - 1;
                }
            }

            for (int i = m; i >= 0; i--)
            {
                if (indexedPeaks[i].ZeroBasedScanIndex < zeroBasedScanIndex)
                {
                    break;
                }

                m--;
            }

            if (m < 0)
            {
                m = 0;
            }

            return m;
        }
        #endregion
    }
}