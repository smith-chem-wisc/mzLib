using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using FlashLFQ.Interfaces;
using Easy.Common.Extensions;

namespace FlashLFQ
{
    public class IndexingEngine<T> : IIndexingEngine where T : IIndexedMzPeak
    {
        internal List<T>[] IndexedPeaks;
        internal const int BinsPerDalton = 100;
        public ScanInfo[] ScanInfoArray { get; private set; }

        /// <summary>
        /// Read in all spectral peaks from the scanArray, index the peaks based on mass and retention time, 
        /// and store them in a jagged array of Lists containing all peaks within a particular mass range
        /// </summary>
        /// <param name="scanArray">An array of raw data scans</param>
        public static IndexingEngine<T> InitializeIndexingEngine(MsDataScan[] scanArray)
        {
            IndexingEngine<T> newEngine = new();
            if(newEngine.IndexPeaks(scanArray))
                return newEngine;
            return null;
        }

        /// <summary>
        /// This factory method returns an IndexingEngine instance where the peaks in all MS1 scans have been indexed. 
        /// This method ignores MS2 scans when indexing
        /// </summary>
        public static IndexingEngine<T> InitializeIndexingEngine(MsDataFile dataFile)
        {
            IndexingEngine<T> newEngine = new();
            var scanArray = dataFile.GetMS1Scans()
                .Where(i => i != null && i.MsnOrder == 1)
                .OrderBy(i => i.OneBasedScanNumber)
                .ToArray();
            if (newEngine.IndexPeaks(scanArray))
                return newEngine;
            return null;
        }

        /// <summary>
        /// Read in all spectral peaks from scans, index the peaks and store them in a list ordered by m/z
        /// </summary>
        /// <param name="scanArray">An array of raw data scans</param>
        /// <param name="scanInfo">Outputs a list of scan information for each scan which is needed for FlashLfq
        public virtual bool IndexPeaks(MsDataScan[] scanArray)
        {
            IndexedPeaks = new List<T>[(int)Math.Ceiling(scanArray.Where(p => p != null
                && p.MassSpectrum.LastX != null).Max(p => p.MassSpectrum.LastX.Value) * BinsPerDalton) + 1];
            ScanInfoArray = new ScanInfo[scanArray.Length];

            for (int scanIndex = 0; scanIndex < scanArray.Length; scanIndex++)
            {
                ScanInfoArray[scanIndex] = new ScanInfo(scanArray[scanIndex].OneBasedScanNumber, scanIndex, scanArray[scanIndex].RetentionTime);

                for (int j = 0; j < scanArray[scanIndex].MassSpectrum.XArray.Length; j++)
                {
                    int roundedMz = (int)Math.Round(scanArray[scanIndex].MassSpectrum.XArray[j] * BinsPerDalton, 0);
                    IndexedPeaks[roundedMz] ??= new List<T>();
                    IndexedPeaks[roundedMz].Add((T)(IIndexedMzPeak)
                        new IndexedMassSpectralPeak(
                            scanArray[scanIndex].MassSpectrum.XArray[j],
                            scanArray[scanIndex].MassSpectrum.YArray[j],
                            scanIndex,
                            scanArray[scanIndex].RetentionTime));
                }
            }
            if (IndexedPeaks == null || IndexedPeaks.Length == 0)
                return false;
            else
                return true;
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
            if (IndexedPeaks == null) throw new MzLibException("Error: Attempt to retrieve indexed peak before peak indexing was performed");
            var bins = GetBinsInRange(mz, ppmTolerance);
            if (bins.Count == 0) return default(T);
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
        /// <param name="retentionTime"> the retention time where peak searching will begin </param>
        /// <param name="missedScansAllowed"> the number of successive missed scans allowed before the xic is terminated </param>
        /// <param name="maxPeakHalfWidth"> the maximum distance from the apex RT of the XIC to both start RT and end RT </param>
        /// <returns> A list of IIndexedMzPeak objects, ordered by retention time </returns>
        public List<IIndexedMzPeak> GetXic(double mz, double retentionTime, PpmTolerance ppmTolerance, int missedScansAllowed, double maxPeakHalfWidth = double.MaxValue)
        {
            // get precursor scan to start at
            int scanIndex = -1;
            foreach (ScanInfo scan in ScanInfoArray)
            {
                if (scan.RetentionTime < retentionTime)
                {
                    scanIndex = scan.ZeroBasedScanIndex;
                }
                else
                {
                    break;
                }
            }

            return GetXic(mz, scanIndex, ppmTolerance, missedScansAllowed, maxPeakHalfWidth);
        }

        /// <summary>
        /// A generic method of peak tracing across the retention time. Finds peaks with a given mz that occur on either side of a given
        /// retention time. Peak searching iterates backwards through the scans until the peak 
        /// is no longer observed (i.e., is absent in more scans than allowed, as defined by the
        /// missedScansAllowed parameter. Missed scans don't have to be sequential. The same procedure
        /// is then repeated in the forward direction.
        /// </summary>
        /// <param name="mz"> the m/z of the peak to be searched for </param>
        /// <param name="zeroBasedStartIndex"> the scan where peak searching begins </param>
        /// <param name="missedScansAllowed"> the number of successive missed scans allowed before the xic is terminated </param>
        /// <param name="maxPeakHalfWidth"> the maximum distance from the apex RT of the XIC to both start RT and end RT </param>
        /// <returns> A list of IIndexedMzPeak objects, ordered by retention time </returns>
        public List<IIndexedMzPeak> GetXic(double mz, int zeroBasedStartIndex, PpmTolerance ppmTolerance, int missedScansAllowed, double maxPeakHalfWidth = double.MaxValue)
        {
            if (IndexedPeaks == null) throw new MzLibException("Error: Attempt to retrieve XIC before peak indexing was performed");
            List<IIndexedMzPeak> xic = new List<IIndexedMzPeak>();
            var allBins = GetBinsInRange(mz, ppmTolerance);
            if (allBins == null || allBins.Count == 0)
                return xic;

            // For each bin, find + store a pointer to the current index
            int[] peakPointerArray = allBins.Select(b => BinarySearchForIndexedPeak(b, zeroBasedStartIndex)).ToArray();
            T initialPeak = GetBestPeakFromBins(allBins, mz, zeroBasedStartIndex, peakPointerArray, ppmTolerance);

            if (initialPeak != null)
                xic.Add(initialPeak);

            foreach (int direction in new List<int> { -1, 1 })
            {
                //int missedPeaks = initialPeak == null ? 1 : 0;
                int missedPeaks = 0; // In the legacy code, the initial peak was not counted as a missed peak if it wasn't found. This should change in the future, but for now, we will keep it as is.
                int currentZeroBasedScanIndex = zeroBasedStartIndex;
                var pointerArrayCopy = new int[peakPointerArray.Length];
                Array.Copy(peakPointerArray, pointerArrayCopy, peakPointerArray.Length);

                while (missedPeaks <= missedScansAllowed)
                {
                    // increment the scan index we're searching for
                    currentZeroBasedScanIndex += direction;
                    if(currentZeroBasedScanIndex < 0 || currentZeroBasedScanIndex > ScanInfoArray.Length - 1)
                        break;
                    
                    for (int i = 0; i < pointerArrayCopy.Length; i++)
                    {
                        // Switch statement designed to increment all pointers in the pointer array
                        // such that they point at the first instances of a peak with the given currentZeroBasedScanIndex
                        switch (direction)
                        {
                            case -1:
                                do
                                {
                                    pointerArrayCopy[i]--;
                                } while (pointerArrayCopy[i] >= 0 && allBins[i][pointerArrayCopy[i]].ZeroBasedScanIndex > currentZeroBasedScanIndex - 1);
                                pointerArrayCopy[i]++;
                                break;
                            case 1:
                                while (pointerArrayCopy[i] < allBins[i].Count - 1 && allBins[i][pointerArrayCopy[i]].ZeroBasedScanIndex < currentZeroBasedScanIndex)
                                    pointerArrayCopy[i]++;
                                break;
                        }
                    }

                    // Search for the next peak
                    T nextPeak = GetBestPeakFromBins(allBins, mz, currentZeroBasedScanIndex, pointerArrayCopy, ppmTolerance);

                    // Add the peak to the XIC or increment the missed peaks
                    if (nextPeak == null)
                        missedPeaks++;
                    else
                    {
                        xic.Add(nextPeak);
                        missedPeaks = 0;
                    }    
                }
            }

            // Sort the XIC in place
            xic.Sort((x, y) => x.RetentionTime.CompareTo(y.RetentionTime));
            return xic;
        }

        #region Peak finding helper functions

        internal List<List<T>> GetBinsInRange(double mz, PpmTolerance ppmTolerance)
        {
            int ceilingMz = (int)Math.Ceiling(ppmTolerance.GetMaximumValue(mz) * BinsPerDalton);
            int floorMz = (int)Math.Floor(ppmTolerance.GetMinimumValue(mz) * BinsPerDalton);
            List<List<T>> allBins = new();
            for (int j = floorMz; j <= ceilingMz; j++)
            {
                if (j >= IndexedPeaks.Length || IndexedPeaks[j] == null)
                    continue;
                allBins.Add(IndexedPeaks[j]);
            }
            return allBins;
        }

        internal static T GetBestPeakFromBins(List<List<T>> allBins, double mz, int zeroBasedScanIndex, IList<int> peakIndicesInBins, PpmTolerance ppmTolerance)
        {
            T bestPeak = default(T);
            for (int i = 0; i < allBins.Count; i++)
            {
                var tempPeak = GetPeakFromBin(allBins[i], mz, zeroBasedScanIndex, peakIndicesInBins[i], ppmTolerance);
                if (tempPeak == null) continue;
                // Check if the peak is within the tolerance and if it is closer to the target Mz than the current peak
                if (bestPeak == null || Math.Abs(tempPeak.Mz - mz) < Math.Abs(bestPeak.Mz - mz))
                {
                    bestPeak = tempPeak;
                }
            }
            return bestPeak;
        }

        /// <summary>
        /// Returns the peak that is closest to the target mz
        /// </summary>
        internal static T GetPeakFromBin(List<T> bin, double mz, int zeroBasedScanIndex, int peakIndexInBin, PpmTolerance ppmTolerance)
        {
            T bestPeak = default(T);
            if (peakIndexInBin < 0 || peakIndexInBin >= bin.Count) return bestPeak;
            for (int i = peakIndexInBin; i < bin.Count; i++)
            {
                T peak = bin[i];

                if (peak.ZeroBasedScanIndex > zeroBasedScanIndex)
                {
                    break;
                }

                if (ppmTolerance.Within(peak.Mz, mz)
                    && peak.ZeroBasedScanIndex == zeroBasedScanIndex
                    && (bestPeak.IsDefault<T>() || Math.Abs(peak.Mz - mz) < Math.Abs(bestPeak.Mz - mz)))
                {
                    bestPeak = peak;
                }
            }
            return bestPeak;
        }

        internal static int BinarySearchForIndexedPeak(List<T> indexedPeaks, int zeroBasedScanIndex)
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
