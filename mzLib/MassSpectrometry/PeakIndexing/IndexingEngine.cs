﻿using Chemistry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

#nullable enable
namespace MassSpectrometry
{
    /// <summary>
    /// IIndexingEngine defines the behaviour needed to efficiently retrieve peaks from a jagged array of indexed peaks.
    /// </summary>
    public abstract class IndexingEngine<T> where T : IIndexedPeak
    {
        /// <summary>
        /// Jagged array. Each index of the array corresponds to a mass bin. Each element of the array is a list of peaks that fall within that mass bin.
        /// Peaks within each mass bin are ordered by scan number, ascending. Due to the width of the mass bin, it is possible to have multiple peaks with the same scan number but different masses in a list
        /// </summary>
        protected List<T>[]? IndexedPeaks;
        protected virtual int BinsPerDalton => 100;
        public ScanInfo[]? ScanInfoArray { get; private set; }

        /// <summary>
        /// Read in all spectral peaks from scans, index the peaks and store them in a list ordered by m/z
        /// </summary>
        /// <param name="scanArray">An array of raw data scans</param>
        public virtual bool IndexPeaks(MsDataScan[] scanArray)
        {
            if(scanArray.IsNullOrEmpty() || scanArray.All(p => p == null))
                return false;

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
                    IndexedPeaks[roundedMz].Add((T)(IIndexedPeak)
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

        /// <summary>
        /// A generic method for finding the closest peak with a specified mass and charge state and in a specified scan. Returns null if no peaks within tolerance are found.
        /// </summary>
        /// <param name="mz"> the m/z of the peak to be searched for </param>
        /// <param name="zeroBasedScanIndex"> the zero based index of the scan where the peak is to be found </param>
        public IIndexedPeak? GetIndexedPeak(double theorMass, int zeroBasedScanIndex, Tolerance ppmTolerance, int chargeState) =>
            GetIndexedPeak(theorMass.ToMz(chargeState), zeroBasedScanIndex, ppmTolerance);

        /// <summary>
        /// A generic method for finding the closest peak with a specified m/z and in a specified scan. Returns null if no peaks within tolerance are found.
        /// </summary>
        /// <param name="mz"> the m/z of the peak to be searched for </param>
        /// <param name="zeroBasedScanIndex"> the zero based index of the scan where the peak is to be found </param>
        public IIndexedPeak? GetIndexedPeak(double mz, int zeroBasedScanIndex, Tolerance ppmTolerance)
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
        /// <returns> A list of IIndexedPeak objects, ordered by retention time </returns>
        public List<IIndexedPeak> GetXic(double mz, double retentionTime, Tolerance ppmTolerance, int missedScansAllowed, double maxPeakHalfWidth = double.MaxValue)
        {
            // get precursor scan to start at
            int scanIndex = -1;
            if (ScanInfoArray == null) throw new MzLibException("Error: Attempt to retrieve XIC before peak indexing was performed");
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
        /// <returns> A list of IIndexedPeak objects, ordered by retention time </returns>
        public List<IIndexedPeak> GetXic(double mz, int zeroBasedStartIndex, Tolerance ppmTolerance, int missedScansAllowed, double maxPeakHalfWidth = double.MaxValue)
        {
            if (IndexedPeaks == null || ScanInfoArray == null) throw new MzLibException("Error: Attempt to retrieve XIC before peak indexing was performed");
            List<IIndexedPeak> xic = new List<IIndexedPeak>();
            var allBins = GetBinsInRange(mz, ppmTolerance);
            if (allBins.Count == 0)
                return xic;

            // For each bin, find + store a pointer to the current index
            int[] peakPointerArray = allBins.Select(b => BinarySearchForIndexedPeak(b, zeroBasedStartIndex)).ToArray();
            var initialPeak = GetBestPeakFromBins(allBins, mz, zeroBasedStartIndex, peakPointerArray, ppmTolerance);

            if (initialPeak.IsNotDefaultOrNull())
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
                    if(currentZeroBasedScanIndex < 0 || currentZeroBasedScanIndex > ScanInfoArray.Length - 1 || (initialPeak.IsNotDefaultOrNull() && Math.Abs(ScanInfoArray[currentZeroBasedScanIndex].RetentionTime - initialPeak.RetentionTime) > maxPeakHalfWidth))
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
                    var nextPeak = GetBestPeakFromBins(allBins, mz, currentZeroBasedScanIndex, pointerArrayCopy, ppmTolerance);

                    // Add the peak to the XIC or increment the missed peaks
                    if (nextPeak == null)
                        missedPeaks++;
                    else
                    {
                        if(initialPeak.IsDefaultOrNull()) initialPeak = nextPeak; // if the initial peak is null, set it to the next peak
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

        internal List<List<T>> GetBinsInRange(double mz, Tolerance ppmTolerance)
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

        internal static T? GetBestPeakFromBins(List<List<T>> allBins, double mz, int zeroBasedScanIndex, IList<int> peakIndicesInBins, Tolerance ppmTolerance)
        {
            T? bestPeak = default(T);
            for (int i = 0; i < allBins.Count; i++)
            {
                var tempPeak = GetPeakFromBin(allBins[i], mz, zeroBasedScanIndex, peakIndicesInBins[i], ppmTolerance);
                if (tempPeak.IsDefaultOrNull()) continue;
                // Check if the peak is within the tolerance and if it is closer to the target M than the current peak
                if (bestPeak.IsDefaultOrNull() || Math.Abs(tempPeak.M - mz) < Math.Abs(bestPeak.M - mz))
                {
                    bestPeak = tempPeak;
                }
            }
            return bestPeak;
        }

        /// <summary>
        /// Returns the peak that is closest to the target mz
        /// </summary>
        internal static T GetPeakFromBin(List<T> bin, double mz, int zeroBasedScanIndex, int peakIndexInBin, Tolerance ppmTolerance)
        {
            T? bestPeak = default(T);
            if (peakIndexInBin < 0 || peakIndexInBin >= bin.Count) return bestPeak;
            for (int i = peakIndexInBin; i < bin.Count; i++)
            {
                var peak = bin[i];

                if (peak.ZeroBasedScanIndex > zeroBasedScanIndex)
                {
                    break;
                }

                if (ppmTolerance.Within(peak.M, mz)
                    && peak.ZeroBasedScanIndex == zeroBasedScanIndex
                    && (bestPeak.IsDefaultOrNull() || Math.Abs(peak.M - mz) < Math.Abs(bestPeak.M - mz)))
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
