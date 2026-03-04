// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Builds a <see cref="DiaScanIndex"/> from an array of <see cref="MsDataScan"/> objects.
    /// 
    /// Phase 16B, Prompt 3: Extended to ingest MS1 scans alongside MS2.
    /// 
    /// MS2 pass (unchanged from Phase 14):
    ///   1. Filter to MsnOrder == 2 scans
    ///   2. Discover isolation windows by clustering IsolationMz / IsolationWidth
    ///   3. Assign integer window IDs
    ///   4. Sort scans within each window by RT
    ///   5. Pack all MS2 peak data into contiguous float arrays
    /// 
    /// MS1 pass (new in Prompt 3):
    ///   1. Filter to MsnOrder == 1 scans
    ///   2. Sort by RetentionTime ascending (required for FindMs1ScanIndexAtRt binary search)
    ///   3. Pack all MS1 peak data into separate contiguous float arrays
    /// 
    /// The builder is intentionally stateless — call Build() as a pure function.
    /// </summary>
    public static class DiaScanIndexBuilder
    {
        // Tolerance for clustering isolation windows by center m/z.
        // Two scans are assigned to the same window if their isolation centers
        // differ by less than this value. Typical DIA windows are 25 m/z wide;
        // center jitter is well under 1 m/z.
        private const double WindowClusteringToleranceMz = 0.5;

        /// <summary>
        /// Builds a <see cref="DiaScanIndex"/> from the provided scans.
        /// 
        /// Both MS1 and MS2 scans are processed in a single call:
        ///   - MS2 scans → indexed by isolation window, sorted by RT within each window
        ///   - MS1 scans → sorted globally by RT, stored in parallel SoA arrays
        /// 
        /// Scans with MsnOrder other than 1 or 2 are silently ignored.
        /// Scans missing IsolationMz or IsolationWidth are silently skipped for MS2
        /// (they cannot be assigned to a window).
        /// </summary>
        /// <param name="scans">All scans from a DIA data file (mix of MS1 and MS2).</param>
        /// <returns>A fully constructed <see cref="DiaScanIndex"/> ready for querying.</returns>
        public static DiaScanIndex Build(MsDataScan[] scans)
        {
            if (scans == null) throw new ArgumentNullException(nameof(scans));

            // ── MS2 pass ────────────────────────────────────────────────────
            var (allMz, allIntensity, scanOffsets, scanLengths, scanWindowIds, scanRts,
                 scanScanNumbers, windowToScanRange, windowLowerBounds, windowUpperBounds)
                = BuildMs2Arrays(scans);

            // ── MS1 pass ────────────────────────────────────────────────────
            var (ms1AllMz, ms1AllIntensity, ms1ScanOffsets, ms1ScanLengths,
                 ms1ScanRts, ms1ScanNumbers)
                = BuildMs1Arrays(scans);

            return new DiaScanIndex(
                allMz, allIntensity,
                scanOffsets, scanLengths, scanWindowIds, scanRts, scanScanNumbers,
                windowToScanRange, windowLowerBounds, windowUpperBounds,
                ms1AllMz, ms1AllIntensity,
                ms1ScanOffsets, ms1ScanLengths, ms1ScanRts, ms1ScanNumbers);
        }

        // ── MS2 builder ─────────────────────────────────────────────────────

        private static (
            float[] allMz,
            float[] allIntensity,
            int[] scanOffsets,
            int[] scanLengths,
            int[] scanWindowIds,
            float[] scanRts,
            int[] scanScanNumbers,
            Dictionary<int, (int Start, int Count)> windowToScanRange,
            float[] windowLowerBounds,
            float[] windowUpperBounds)
        BuildMs2Arrays(MsDataScan[] allScans)
        {
            // Filter to valid MS2 scans (must have isolation window metadata)
            var ms2Scans = new List<MsDataScan>(allScans.Length);
            foreach (var scan in allScans)
            {
                if (scan.MsnOrder == 2 && scan.IsolationMz.HasValue && scan.IsolationWidth.HasValue)
                    ms2Scans.Add(scan);
            }

            // ── Step 1: Discover isolation windows ──────────────────────────
            // Two scans belong to the same window if their isolation centers are within
            // WindowClusteringToleranceMz of each other. We sort by isolation center and
            // then do a single linear pass to group them.
            //
            // Window ID is the index in the sorted unique center list.
            var isolationCenters = ms2Scans
                .Select(s => s.IsolationMz!.Value)
                .Distinct()
                .OrderBy(c => c)
                .ToList();

            // Cluster nearby centers (handles minor instrument jitter)
            var clusteredCenters = new List<double>();
            foreach (double center in isolationCenters)
            {
                if (clusteredCenters.Count == 0 ||
                    center - clusteredCenters[clusteredCenters.Count - 1] > WindowClusteringToleranceMz)
                {
                    clusteredCenters.Add(center);
                }
                // else: merge into previous cluster (represented by its first member)
            }

            int windowCount = clusteredCenters.Count;

            // Build window bound arrays: center ± (width / 2)
            // We derive the representative width for each window from the first scan assigned to it.
            var windowWidthMap = new Dictionary<int, double>(windowCount);
            var windowLowerBounds = new float[windowCount];
            var windowUpperBounds = new float[windowCount];

            // Assign each scan a window ID
            int[] scanWindowAssignments = new int[ms2Scans.Count];
            for (int i = 0; i < ms2Scans.Count; i++)
            {
                double center = ms2Scans[i].IsolationMz!.Value;
                int wid = FindClusterIndex(clusteredCenters, center, WindowClusteringToleranceMz);
                scanWindowAssignments[i] = wid;

                // Record window width from first scan in this window
                if (!windowWidthMap.ContainsKey(wid))
                {
                    double width = ms2Scans[i].IsolationWidth!.Value;
                    windowWidthMap[wid] = width;
                    windowLowerBounds[wid] = (float)(clusteredCenters[wid] - width / 2.0);
                    windowUpperBounds[wid] = (float)(clusteredCenters[wid] + width / 2.0);
                }
            }

            // ── Step 2: Group scans by window, sort by RT within each window ─
            // Build per-window lists of (RT, scan list index)
            var windowScanLists = new List<(float Rt, int ScanListIdx)>[windowCount];
            for (int w = 0; w < windowCount; w++)
                windowScanLists[w] = new List<(float, int)>();

            for (int i = 0; i < ms2Scans.Count; i++)
                windowScanLists[scanWindowAssignments[i]].Add(((float)ms2Scans[i].RetentionTime, i));

            foreach (var list in windowScanLists)
                list.Sort((a, b) => a.Rt.CompareTo(b.Rt));

            // ── Step 3: Build the ordered scan list and SoA arrays ──────────
            // Total peaks across all MS2 scans
            long totalPeaks = 0;
            foreach (var scan in ms2Scans)
                totalPeaks += scan.MassSpectrum?.Size ?? 0;

            if (totalPeaks > int.MaxValue)
                throw new InvalidOperationException(
                    $"Total MS2 peak count ({totalPeaks:N0}) exceeds int.MaxValue. " +
                    "Split the file into smaller segments.");

            var allMz = new float[(int)totalPeaks];
            var allIntensity = new float[(int)totalPeaks];
            var scanOffsets = new int[ms2Scans.Count];
            var scanLengths = new int[ms2Scans.Count];
            var scanWindowIds = new int[ms2Scans.Count];
            var scanRts = new float[ms2Scans.Count];
            var scanScanNumbers = new int[ms2Scans.Count];
            var windowToScanRange = new Dictionary<int, (int Start, int Count)>(windowCount);

            int writeOffset = 0;
            int globalScanIdx = 0;

            for (int w = 0; w < windowCount; w++)
            {
                var list = windowScanLists[w];
                int windowStart = globalScanIdx;

                foreach (var (rt, srcIdx) in list)
                {
                    var scan = ms2Scans[srcIdx];
                    var spectrum = scan.MassSpectrum;
                    int peakCount = spectrum?.Size ?? 0;

                    scanOffsets[globalScanIdx] = writeOffset;
                    scanLengths[globalScanIdx] = peakCount;
                    scanWindowIds[globalScanIdx] = w;
                    scanRts[globalScanIdx] = (float)scan.RetentionTime;
                    scanScanNumbers[globalScanIdx] = scan.OneBasedScanNumber;

                    if (peakCount > 0)
                    {
                        double[] xArray = spectrum!.XArray;
                        double[] yArray = spectrum.YArray;
                        for (int p = 0; p < peakCount; p++)
                        {
                            allMz[writeOffset + p] = (float)xArray[p];
                            allIntensity[writeOffset + p] = (float)yArray[p];
                        }
                        writeOffset += peakCount;
                    }

                    globalScanIdx++;
                }

                windowToScanRange[w] = (windowStart, list.Count);
            }

            return (allMz, allIntensity, scanOffsets, scanLengths, scanWindowIds, scanRts,
                    scanScanNumbers, windowToScanRange, windowLowerBounds, windowUpperBounds);
        }

        // ── MS1 builder ─────────────────────────────────────────────────────

        /// <summary>
        /// Builds the MS1 SoA arrays from the full scan array.
        /// 
        /// Steps:
        ///   1. Filter to MsnOrder == 1
        ///   2. Sort by RetentionTime ascending
        ///   3. Pack all peak data into contiguous float arrays
        /// 
        /// No window assignment — MS1 scans are not associated with isolation windows.
        /// </summary>
        private static (
            float[] ms1AllMz,
            float[] ms1AllIntensity,
            int[] ms1ScanOffsets,
            int[] ms1ScanLengths,
            float[] ms1ScanRts,
            int[] ms1ScanNumbers)
        BuildMs1Arrays(MsDataScan[] allScans)
        {
            // Filter and sort by RT
            var ms1Scans = allScans
                .Where(s => s.MsnOrder == 1)
                .OrderBy(s => s.RetentionTime)
                .ToArray();

            int scanCount = ms1Scans.Length;

            if (scanCount == 0)
            {
                // Return empty but non-null arrays so DiaScanIndex constructor
                // receives valid inputs and Ms1ScanCount == 0.
                return (Array.Empty<float>(), Array.Empty<float>(),
                        Array.Empty<int>(), Array.Empty<int>(),
                        Array.Empty<float>(), Array.Empty<int>());
            }

            // Count total MS1 peaks
            long totalPeaks = 0;
            foreach (var scan in ms1Scans)
                totalPeaks += scan.MassSpectrum?.Size ?? 0;

            if (totalPeaks > int.MaxValue)
                throw new InvalidOperationException(
                    $"Total MS1 peak count ({totalPeaks:N0}) exceeds int.MaxValue.");

            var ms1AllMz = new float[(int)totalPeaks];
            var ms1AllIntensity = new float[(int)totalPeaks];
            var ms1ScanOffsets = new int[scanCount];
            var ms1ScanLengths = new int[scanCount];
            var ms1ScanRts = new float[scanCount];
            var ms1ScanNumbers = new int[scanCount];

            int writeOffset = 0;

            for (int i = 0; i < scanCount; i++)
            {
                var scan = ms1Scans[i];
                var spectrum = scan.MassSpectrum;
                int peakCount = spectrum?.Size ?? 0;

                ms1ScanOffsets[i] = writeOffset;
                ms1ScanLengths[i] = peakCount;
                ms1ScanRts[i] = (float)scan.RetentionTime;
                ms1ScanNumbers[i] = scan.OneBasedScanNumber;

                if (peakCount > 0)
                {
                    double[] xArray = spectrum!.XArray;
                    double[] yArray = spectrum.YArray;
                    for (int p = 0; p < peakCount; p++)
                    {
                        ms1AllMz[writeOffset + p] = (float)xArray[p];
                        ms1AllIntensity[writeOffset + p] = (float)yArray[p];
                    }
                    writeOffset += peakCount;
                }
            }

            return (ms1AllMz, ms1AllIntensity, ms1ScanOffsets, ms1ScanLengths,
                    ms1ScanRts, ms1ScanNumbers);
        }

        // ── Utilities ────────────────────────────────────────────────────────

        /// <summary>
        /// Binary search for the cluster index in the sorted clusteredCenters list
        /// such that |clusteredCenters[idx] - center| is minimized.
        /// Falls back to linear scan for the rare case of exact ties at boundaries.
        /// </summary>
        private static int FindClusterIndex(
            List<double> clusteredCenters, double center, double tolerance)
        {
            int lo = 0, hi = clusteredCenters.Count - 1;
            while (lo < hi)
            {
                int mid = lo + ((hi - lo) >> 1);
                if (clusteredCenters[mid] < center - tolerance)
                    lo = mid + 1;
                else
                    hi = mid;
            }
            // lo is the first cluster center within tolerance
            // If it's still more than tolerance away, snap to nearest
            if (lo > 0 && Math.Abs(center - clusteredCenters[lo - 1]) <
                          Math.Abs(center - clusteredCenters[lo]))
                return lo - 1;
            return lo;
        }
    }
}