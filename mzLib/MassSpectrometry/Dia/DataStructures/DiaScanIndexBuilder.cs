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
    /// Window Discovery:
    ///   Real DIA instruments (Thermo, Bruker, SCIEX) report isolation window boundaries
    ///   with slight m/z jitter (±0.1 Da scan-to-scan) due to calibration drift, rounding,
    ///   or firmware differences. A typical DIA scheme has 18-40 windows spanning 25-110 Da
    ///   each, so the jitter is far smaller than the gap between true windows.
    ///   
    ///   We use a Da-based tolerance (default 1.0 Da) on the isolation window center to
    ///   group scans into windows. This absorbs all realistic instrument jitter while never
    ///   merging truly distinct windows (which differ by 25+ Da at minimum in any DIA scheme).
    /// 
    /// Usage:
    ///   var scans = myMsDataFile.GetMsDataScans();
    ///   var index = DiaScanIndexBuilder.Build(scans);
    /// </summary>
    public static class DiaScanIndexBuilder
    {
        /// <summary>
        /// Default Da tolerance for grouping isolation window centers.
        /// Real DIA windows differ by 25+ Da; instrument jitter is &lt;0.5 Da.
        /// 1.0 Da safely absorbs all jitter without merging distinct windows.
        /// </summary>
        public const double DefaultWindowGroupingToleranceDa = 1.0;

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
        /// <param name="isolationWindowGroupingToleranceDa">
        /// Dalton tolerance for grouping isolation window centers into the same window.
        /// Default 1.0 Da. Two scans whose isolation centers differ by less than this 
        /// tolerance are assigned the same window ID. This should be large enough to
        /// absorb instrument jitter (typically ±0.1 Da) but small enough to never merge
        /// distinct DIA windows (typically 25+ Da apart).
        /// </param>
        /// <returns>A fully built DiaScanIndex ready for fragment extraction.</returns>
        public static DiaScanIndex Build(MsDataScan[] allScans, double isolationWindowGroupingToleranceDa = DefaultWindowGroupingToleranceDa)
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
            var windowDefinitions = DiscoverIsolationWindows(ms2Scans, isolationWindowGroupingToleranceDa);
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
        /// Overload that accepts the legacy PPM parameter name for backward compatibility
        /// with existing tests. Converts to a reasonable Da tolerance internally.
        /// 
        /// For new code, prefer the default Build(allScans) or Build(allScans, toleranceDa).
        /// </summary>
        [Obsolete("Use Build(allScans, isolationWindowGroupingToleranceDa) instead. PPM-based grouping is too tight for real instrument data.")]
        public static DiaScanIndex BuildWithPpmTolerance(MsDataScan[] allScans, double isolationWindowTolerancePpm)
        {
            // Convert PPM at a typical DIA center (~600 m/z) to Da, with a floor of 0.5 Da
            double daTolerance = Math.Max(0.5, 600.0 * isolationWindowTolerancePpm / 1e6);
            return Build(allScans, daTolerance);
        }

        /// <summary>
        /// Discovers unique DIA isolation windows from a set of MS2 scans.
        /// 
        /// Algorithm:
        ///   1. Collect all (center, width) pairs from scans
        ///   2. Sort by center m/z ascending
        ///   3. Greedy merge: walk sorted centers, start a new cluster whenever the center
        ///      is more than toleranceDa away from the running cluster mean
        ///   4. For each cluster, compute the consensus window bounds from the most common
        ///      isolation width (mode), applied symmetrically around the mean center
        /// 
        /// This approach handles:
        ///   - Instrument jitter (±0.1 Da typical on Thermo QE/Exploris)
        ///   - Slight rounding differences in mzML conversion
        ///   - Missing IsolationWidth on some scans (uses cluster consensus)
        /// 
        /// Returns windows sorted by lower bound (ascending).
        /// </summary>
        private static List<WindowDefinition> DiscoverIsolationWindows(
            List<MsDataScan> ms2Scans, double toleranceDa)
        {
            // ── Collect all (center, width) observations ────────────────────
            // We store the raw values for every scan so we can compute robust statistics
            var observations = new List<(double Center, double Width)>(ms2Scans.Count);
            for (int i = 0; i < ms2Scans.Count; i++)
            {
                double center = ms2Scans[i].IsolationMz.Value;
                double width = ms2Scans[i].IsolationWidth.HasValue
                    ? ms2Scans[i].IsolationWidth.Value
                    : 0.0;
                observations.Add((center, width));
            }

            // Sort by center ascending for greedy clustering
            observations.Sort((a, b) => a.Center.CompareTo(b.Center));

            // ── Greedy clustering ───────────────────────────────────────────
            // Walk sorted observations, group consecutive scans whose centers are
            // within toleranceDa of the running cluster mean.
            var clusters = new List<List<(double Center, double Width)>>();
            var currentCluster = new List<(double Center, double Width)> { observations[0] };
            double clusterSum = observations[0].Center;

            for (int i = 1; i < observations.Count; i++)
            {
                double clusterMean = clusterSum / currentCluster.Count;
                if (Math.Abs(observations[i].Center - clusterMean) <= toleranceDa)
                {
                    // Same window cluster
                    currentCluster.Add(observations[i]);
                    clusterSum += observations[i].Center;
                }
                else
                {
                    // New cluster — save current and start fresh
                    clusters.Add(currentCluster);
                    currentCluster = new List<(double Center, double Width)> { observations[i] };
                    clusterSum = observations[i].Center;
                }
            }
            clusters.Add(currentCluster); // Don't forget the last cluster

            // ── Build window definitions from clusters ──────────────────────
            var windows = new List<WindowDefinition>(clusters.Count);
            for (int c = 0; c < clusters.Count; c++)
            {
                var cluster = clusters[c];

                // Consensus center = mean of all observations in cluster
                double centerSum = 0;
                for (int i = 0; i < cluster.Count; i++)
                    centerSum += cluster[i].Center;
                double consensusCenter = centerSum / cluster.Count;

                // Consensus width = most common non-zero width (mode), or fallback to max
                double consensusWidth = ComputeConsensusWidth(cluster);

                double halfWidth = consensusWidth / 2.0;
                windows.Add(new WindowDefinition(
                    (float)(consensusCenter - halfWidth),
                    (float)(consensusCenter + halfWidth),
                    (float)consensusCenter));
            }

            // Sort by lower bound (should already be sorted, but be safe)
            windows.Sort((a, b) => a.LowerBound.CompareTo(b.LowerBound));

            return windows;
        }

        /// <summary>
        /// Computes the consensus isolation width from a cluster of observations.
        /// Uses the mode (most frequent) non-zero width, rounded to 0.1 Da to absorb jitter.
        /// Falls back to the maximum width if no non-zero widths are present.
        /// </summary>
        private static double ComputeConsensusWidth(List<(double Center, double Width)> cluster)
        {
            // Round widths to 0.1 Da and count occurrences
            var widthCounts = new Dictionary<double, int>();
            double maxWidth = 0;

            for (int i = 0; i < cluster.Count; i++)
            {
                double w = cluster[i].Width;
                if (w <= 0) continue;
                if (w > maxWidth) maxWidth = w;

                // Round to nearest 0.1 to group near-identical widths
                double rounded = Math.Round(w, 1);
                if (widthCounts.ContainsKey(rounded))
                    widthCounts[rounded]++;
                else
                    widthCounts[rounded] = 1;
            }

            if (widthCounts.Count == 0)
                return 0.0; // No width information available

            // Return the most common rounded width
            double bestWidth = 0;
            int bestCount = 0;
            foreach (var kvp in widthCounts)
            {
                if (kvp.Value > bestCount)
                {
                    bestCount = kvp.Value;
                    bestWidth = kvp.Key;
                }
            }

            return bestWidth;
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