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
        private List<IndexedMassSpectralPeak>[] _indexedPeaks;
        private readonly Serializer _serializer;
        private const int BinsPerDalton = 100;
        public Ms1ScanInfo[] ScanInfoArray { get; private set; }
        public SpectraFileInfo SpectraFile { get; init; }

        public PeakIndexingEngine(SpectraFileInfo file)
        {
            var messageTypes = new List<Type>
            {
                typeof(List<IndexedMassSpectralPeak>[]), typeof(List<IndexedMassSpectralPeak>),
                typeof(IndexedMassSpectralPeak)
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
            _indexedPeaks = new List<IndexedMassSpectralPeak>[(int)Math.Ceiling(msDataScans.Where(p => p != null
                && p.MassSpectrum.LastX != null).Max(p => p.MassSpectrum.LastX.Value) * BinsPerDalton) + 1];
            ScanInfoArray = new Ms1ScanInfo[msDataScans.Length];

            for (int scanIndex = 0; scanIndex < msDataScans.Length; scanIndex++)
            {
                ScanInfoArray[scanIndex] = new Ms1ScanInfo(msDataScans[scanIndex].OneBasedScanNumber, scanIndex, msDataScans[scanIndex].RetentionTime);

                for (int j = 0; j < msDataScans[scanIndex].MassSpectrum.XArray.Length; j++)
                {
                    int roundedMz = (int)Math.Round(msDataScans[scanIndex].MassSpectrum.XArray[j] * BinsPerDalton, 0);
                    _indexedPeaks[roundedMz] ??= new List<IndexedMassSpectralPeak>();
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
                _indexedPeaks = (List<IndexedMassSpectralPeak>[])_serializer.Deserialize(indexFile);
            }

            File.Delete(indexPath);
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
        public  List<IIndexedMzPeak> GetXIC(double mz, int zeroBasedStartIndex, PpmTolerance ppmTolerance, int missedScansAllowed, double maxPeakHalfWidth = double.MaxValue)
        {
            var xic = new List<IIndexedMzPeak>();

            // go right
            int missedScans = 0;
            for (int t = zeroBasedStartIndex; t < ScanInfoArray.Length; t++)
            {
                var peak = GetIndexedPeak(mz, t, ppmTolerance);

                if (peak == null && t != zeroBasedStartIndex)
                {
                    missedScans++;
                }
                else if (peak != null)
                {
                    missedScans = 0;
                    xic.Add(peak);

                    if (peak.RetentionTime - xic.First().RetentionTime > maxPeakHalfWidth)
                    {
                        break;
                    }
                }

                if (missedScans > missedScansAllowed)
                {
                    break;
                }
            }

            // go left
            missedScans = 0;
            for (int t = zeroBasedStartIndex - 1; t >= 0; t--)
            {
                var peak = GetIndexedPeak(mz, t, ppmTolerance);

                if (peak == null && t != zeroBasedStartIndex)
                {
                    missedScans++;
                }
                else if (peak != null)
                {
                    missedScans = 0;
                    xic.Add(peak);

                    if (xic.First().RetentionTime - peak.RetentionTime > maxPeakHalfWidth)
                    {
                        break;
                    }
                }

                if (missedScans > missedScansAllowed)
                {
                    break;
                }
            }

            // Sorts the list by RT in-place. (OrderBy does not sort in place, but creates a copy instead)
            xic.Sort((x, y) => x.RetentionTime.CompareTo(y.RetentionTime));

            return xic;
        }


        public IIndexedMzPeak GetIndexedPeak(double theorMass, int zeroBasedScanIndex, PpmTolerance ppmTolerance, int chargeState) =>
            GetIndexedPeak(theorMass.ToMz(chargeState), zeroBasedScanIndex, ppmTolerance);

        /// <summary>
        /// A generic method for finding the closest peak with a specified m/z and in a specified scan. Returns null if no peaks within tolerance are found.
        /// </summary>
        public IIndexedMzPeak GetIndexedPeak(double mz, int zeroBasedScanIndex, PpmTolerance ppmTolerance)
        {
            IndexedMassSpectralPeak bestPeak = null;
            int ceilingMz = (int)Math.Ceiling(ppmTolerance.GetMaximumValue(mz) * BinsPerDalton);
            int floorMz = (int)Math.Floor(ppmTolerance.GetMinimumValue(mz) * BinsPerDalton);

            for (int j = floorMz; j <= ceilingMz; j++)
            {
                if (j < _indexedPeaks.Length && _indexedPeaks[j] != null)
                {
                    List<IndexedMassSpectralPeak> bin = _indexedPeaks[j];
                    int index = BinarySearchForIndexedPeak(bin, zeroBasedScanIndex);

                    for (int i = index; i < bin.Count; i++)
                    {
                        IndexedMassSpectralPeak peak = bin[i];

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
                }
            }

            return bestPeak;
        }

        public IIndexedMzPeak GetPeakFromBin(double mz, int zeroBasedScanIndex, int binIndex, int peakListIndex, PpmTolerance ppmTolerance)
        {
            
            if (binIndex >= _indexedPeaks.Length || _indexedPeaks[binIndex] == null)
                return null;
            List<IndexedMassSpectralPeak> bin = _indexedPeaks[binIndex];

            IndexedMassSpectralPeak bestPeak = null;
            for (int i = peakListIndex; i < bin.Count; i++)
            {
                IndexedMassSpectralPeak peak = bin[i];

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

        public List<IIndexedMzPeak> GetXicFromBin(double mz, int zeroBasedScanIndex, int binIndex, PpmTolerance ppmTolerance)
        {
            List<IIndexedMzPeak> xic = new List<IIndexedMzPeak>();
            IndexedMassSpectralPeak bestPeak = null;
            if (binIndex >= _indexedPeaks.Length || _indexedPeaks[binIndex] == null)
                return null;

            List<IndexedMassSpectralPeak> bin = _indexedPeaks[binIndex];
            int peakListIndex = BinarySearchForIndexedPeak(bin, zeroBasedScanIndex);
            IIndexedMzPeak initialPeak = GetPeakFromBin(mz, zeroBasedScanIndex, binIndex, peakListIndex, ppmTolerance);
            if ( initialPeak == null)
            {
                return xic;
            }

            foreach (int direction in new List<int> { -1, 1})
            {
                for (int i = peakListIndex + direction; i < bin.Count && i >= 0; i+=direction)
                {
                    IIndexedMzPeak nextPeak = GetPeakFromBin(mz, zeroBasedScanIndex + direction, binIndex, i, ppmTolerance);
                    IndexedMassSpectralPeak peak = bin[i];

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
            }
            
            
            return bestPeak;
        }

        private static int BinarySearchForIndexedPeak(List<IndexedMassSpectralPeak> indexedPeaks, int zeroBasedScanIndex)
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
    }
}