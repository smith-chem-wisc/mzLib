// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Builds a DiaScanIndex from an array of MsDataScan objects (as loaded by mzLib's
    /// existing Mzml or ThermoRawFileReader classes).
    /// 
    /// Design rationale:
    ///   - This is the bridge between mzLib's object-oriented scan model and the 
    ///     performance-oriented SoA layout required for DIA processing.
    ///   - The builder iterates the scan array twice: once to count total peaks (so we can
    ///     allocate exact-sized arrays), and once to copy peak data into the flat layout.
    ///   - After building, the original MsDataScan objects can be released by the caller.
    ///   - Window IDs are assigned by sorting unique isolation window boundaries and mapping
    ///     each scan's isolation center to the nearest window definition.
    /// 
    /// Usage:
    ///   var scans = myMsDataFile.GetMsDataScans();
    ///   var index = DiaScanIndexBuilder.Build(scans);
    /// </summary>
    public static class DiaScanIndexBuilder
    {
        /// <summary>
        /// Builds a DiaScanIndex from an array of MsDataScan objects.
        /// Only MS2 (MsnOrder == 2) centroided scans with valid isolation windows are included.
        /// 
        /// Steps:
        ///   1. Filter to MS2 scans with valid isolation info
        ///   2. Discover unique isolation windows and assign integer IDs
        ///   3. Sort scans by (windowId, RT) so scans within each window are RT-ordered
        ///   4. Count total peaks across all qualifying scans
        ///   5. Allocate contiguous float arrays and copy peak data (double → float)
        ///   6. Build the window → scan range index
        /// </summary>
        /// <param name="allScans">All scans from an MsDataFile (MS1 + MS2 mixed).</param>
        /// <param name="isolationWindowTolerancePpm">
        /// PPM tolerance for grouping isolation window centers into the same window.
        /// Default 10 ppm. Two scans whose isolation centers differ by less than this 
        /// tolerance are assigned the same window ID.
        /// </param>
        /// <returns>A fully built DiaScanIndex ready for fragment extraction.</returns>
        public static DiaScanIndex Build(MsDataScan[] allScans, double isolationWindowTolerancePpm = 10.0)
        {
            if (allScans == null) throw new ArgumentNullException(nameof(allScans));

            // ── Step 1: Filter to valid MS2 scans ───────────────────────────
            var ms2Scans = new List<MsDataScan>();
            for (int i = 0; i < allScans.Length; i++)
            {
                var scan = allScans[i];
                if (scan == null) continue;
                if (scan.MsnOrder != 2) continue;
                if (!scan.IsolationMz.HasValue) continue;
                if (scan.MassSpectrum == null || scan.MassSpectrum.Size == 0) continue;

                ms2Scans.Add(scan);
            }

            if (ms2Scans.Count == 0)
            {
                return BuildEmpty();
            }

            // ── Step 2: Discover unique windows and assign IDs ──────────────
            var windowDefinitions = DiscoverIsolationWindows(ms2Scans, isolationWindowTolerancePpm);
            // windowDefinitions: sorted list of (lowerBound, upperBound, centroidMz)
            // Window ID = index in this list

            // Map each scan to its window ID
            var scanWindowIds = new int[ms2Scans.Count];
            for (int i = 0; i < ms2Scans.Count; i++)
            {
                scanWindowIds[i] = AssignWindowId(ms2Scans[i].IsolationMz.Value, windowDefinitions);
            }

            // ── Step 3: Sort scans by (windowId, RT) ────────────────────────
            // Create index array and sort it so we don't move the actual scan objects
            var sortedIndices = Enumerable.Range(0, ms2Scans.Count).ToArray();
            Array.Sort(sortedIndices, (a, b) =>
            {
                int windowCmp = scanWindowIds[a].CompareTo(scanWindowIds[b]);
                if (windowCmp != 0) return windowCmp;
                return ms2Scans[a].RetentionTime.CompareTo(ms2Scans[b].RetentionTime);
            });

            // ── Step 4: Count total peaks ───────────────────────────────────
            long totalPeaks = 0;
            for (int i = 0; i < ms2Scans.Count; i++)
            {
                totalPeaks += ms2Scans[i].MassSpectrum.Size;
            }

            // Guard against exceeding int.MaxValue (would need chunked layout for >2B peaks)
            if (totalPeaks > int.MaxValue)
                throw new InvalidOperationException(
                    $"Total peak count ({totalPeaks}) exceeds int.MaxValue. File is too large for single-index layout.");

            // ── Step 5: Allocate and fill contiguous arrays ─────────────────
            int totalPeakCount = (int)totalPeaks;
            float[] allMz = new float[totalPeakCount];
            float[] allIntensity = new float[totalPeakCount];

            int scanCount = ms2Scans.Count;
            int[] offsets = new int[scanCount];
            int[] lengths = new int[scanCount];
            int[] windowIds = new int[scanCount];
            float[] rts = new float[scanCount];
            int[] oneBasedScanNumbers = new int[scanCount];

            int currentOffset = 0;
            for (int outputIdx = 0; outputIdx < scanCount; outputIdx++)
            {
                int sourceIdx = sortedIndices[outputIdx];
                var scan = ms2Scans[sourceIdx];
                var spectrum = scan.MassSpectrum;
                int peakCount = spectrum.Size;

                offsets[outputIdx] = currentOffset;
                lengths[outputIdx] = peakCount;
                windowIds[outputIdx] = scanWindowIds[sourceIdx];
                rts[outputIdx] = (float)scan.RetentionTime;
                oneBasedScanNumbers[outputIdx] = scan.OneBasedScanNumber;

                // Copy peaks: double → float conversion
                double[] mzSource = spectrum.XArray;
                double[] intSource = spectrum.YArray;
                for (int p = 0; p < peakCount; p++)
                {
                    allMz[currentOffset + p] = (float)mzSource[p];
                    allIntensity[currentOffset + p] = (float)intSource[p];
                }

                currentOffset += peakCount;
            }

            // ── Step 6: Build window → scan range index ─────────────────────
            var windowToScanRange = new Dictionary<int, (int Start, int Count)>();
            if (scanCount > 0)
            {
                int rangeStart = 0;
                int currentWindowId = windowIds[0];
                for (int i = 1; i < scanCount; i++)
                {
                    if (windowIds[i] != currentWindowId)
                    {
                        windowToScanRange[currentWindowId] = (rangeStart, i - rangeStart);
                        currentWindowId = windowIds[i];
                        rangeStart = i;
                    }
                }
                // Final window
                windowToScanRange[currentWindowId] = (rangeStart, scanCount - rangeStart);
            }

            // Build window bound arrays
            float[] windowLower = new float[windowDefinitions.Count];
            float[] windowUpper = new float[windowDefinitions.Count];
            for (int i = 0; i < windowDefinitions.Count; i++)
            {
                windowLower[i] = windowDefinitions[i].LowerBound;
                windowUpper[i] = windowDefinitions[i].UpperBound;
            }

            return new DiaScanIndex(
                allMz, allIntensity,
                offsets, lengths, windowIds, rts, oneBasedScanNumbers,
                windowToScanRange,
                windowLower, windowUpper);
        }

        /// <summary>
        /// Discovers unique DIA isolation windows from a set of MS2 scans.
        /// 
        /// Groups scans by their isolation center m/z (within the given ppm tolerance),
        /// then computes the window bounds from IsolationMz ± IsolationWidth/2.
        /// If IsolationWidth is not available, uses just the center m/z as a point window.
        /// 
        /// Returns windows sorted by lower bound (ascending).
        /// </summary>
        private static List<WindowDefinition> DiscoverIsolationWindows(
            List<MsDataScan> ms2Scans, double tolerancePpm)
        {
            // Collect all unique isolation centers
            var centers = new List<double>();
            var widths = new Dictionary<double, double>(); // center → width (from first scan that has it)

            for (int i = 0; i < ms2Scans.Count; i++)
            {
                double center = ms2Scans[i].IsolationMz.Value;
                bool isNew = true;

                for (int j = 0; j < centers.Count; j++)
                {
                    double ppmDiff = Math.Abs(center - centers[j]) / centers[j] * 1e6;
                    if (ppmDiff < tolerancePpm)
                    {
                        isNew = false;
                        break;
                    }
                }

                if (isNew)
                {
                    centers.Add(center);
                    if (ms2Scans[i].IsolationWidth.HasValue)
                    {
                        widths[center] = ms2Scans[i].IsolationWidth.Value;
                    }
                }
            }

            // Sort centers ascending
            centers.Sort();

            // Build window definitions
            var windows = new List<WindowDefinition>(centers.Count);
            for (int i = 0; i < centers.Count; i++)
            {
                double center = centers[i];
                double halfWidth = widths.ContainsKey(center) ? widths[center] / 2.0 : 0.0;
                windows.Add(new WindowDefinition(
                    (float)(center - halfWidth),
                    (float)(center + halfWidth),
                    (float)center));
            }

            return windows;
        }

        /// <summary>
        /// Assigns a window ID (index into the sorted window definitions list) to a scan
        /// based on its isolation center m/z. Finds the nearest window center.
        /// </summary>
        private static int AssignWindowId(double isolationMz, List<WindowDefinition> windows)
        {
            int bestId = 0;
            double bestDist = Math.Abs(isolationMz - windows[0].CenterMz);

            for (int i = 1; i < windows.Count; i++)
            {
                double dist = Math.Abs(isolationMz - windows[i].CenterMz);
                if (dist < bestDist)
                {
                    bestDist = dist;
                    bestId = i;
                }
            }

            return bestId;
        }

        /// <summary>Returns an empty DiaScanIndex (no scans, no windows).</summary>
        private static DiaScanIndex BuildEmpty()
        {
            return new DiaScanIndex(
                Array.Empty<float>(), Array.Empty<float>(),
                Array.Empty<int>(), Array.Empty<int>(), Array.Empty<int>(),
                Array.Empty<float>(), Array.Empty<int>(),
                new Dictionary<int, (int, int)>(),
                Array.Empty<float>(), Array.Empty<float>());
        }

        /// <summary>
        /// Represents a discovered DIA isolation window with its m/z boundaries and center.
        /// </summary>
        internal readonly struct WindowDefinition
        {
            public readonly float LowerBound;
            public readonly float UpperBound;
            public readonly float CenterMz;

            public WindowDefinition(float lowerBound, float upperBound, float centerMz)
            {
                LowerBound = lowerBound;
                UpperBound = upperBound;
                CenterMz = centerMz;
            }
        }
    }
}
