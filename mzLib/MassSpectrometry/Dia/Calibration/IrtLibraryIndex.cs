// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Collections.Generic;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Sorted index of library precursors by iRT value for O(log N + K) candidate selection.
    /// 
    /// The index sorts LibraryPrecursorInput entries by their IrtValue and stores parallel arrays
    /// for binary search. Given a calibrated iRT and window half-width, QueryByIrt returns
    /// the subset of precursors whose iRT falls within the window.
    /// 
    /// This is the core data structure for the "refined candidate selection" step (Section 5
    /// of the iRT calibration spec). It avoids the previous approach of applying a fixed
    /// RT tolerance per-precursor and instead enables scan-centric querying: for each scan,
    /// convert RT to iRT, then binary search the library for candidates.
    /// 
    /// Performance: O(log N) for the binary search + O(K) for iterating candidates,
    /// where N = total library size and K = candidates in the window.
    /// 
    /// Immutable after construction. Thread-safe for concurrent reads.
    /// </summary>
    public sealed class IrtLibraryIndex
    {
        // Sorted iRT values (parallel to _sortedIndices)
        private readonly double[] _sortedIrts;

        // Original indices into the input precursor list, sorted by iRT
        private readonly int[] _sortedIndices;

        /// <summary>Number of precursors in the index.</summary>
        public int Count => _sortedIrts.Length;

        /// <summary>Minimum iRT in the library. NaN if empty.</summary>
        public double MinIrt => Count > 0 ? _sortedIrts[0] : double.NaN;

        /// <summary>Maximum iRT in the library. NaN if empty.</summary>
        public double MaxIrt => Count > 0 ? _sortedIrts[Count - 1] : double.NaN;

        /// <summary>
        /// Builds a sorted iRT index from a list of precursors.
        /// Only precursors with non-null IrtValue are included.
        /// </summary>
        /// <param name="precursors">The library precursor list.</param>
        public IrtLibraryIndex(IList<LibraryPrecursorInput> precursors)
        {
            if (precursors == null) throw new ArgumentNullException(nameof(precursors));

            // Count precursors with valid iRT
            int countWithIrt = 0;
            for (int i = 0; i < precursors.Count; i++)
            {
                if (precursors[i].IrtValue.HasValue)
                    countWithIrt++;
            }

            _sortedIrts = new double[countWithIrt];
            _sortedIndices = new int[countWithIrt];

            // Fill
            int idx = 0;
            for (int i = 0; i < precursors.Count; i++)
            {
                if (precursors[i].IrtValue.HasValue)
                {
                    _sortedIrts[idx] = precursors[i].IrtValue.Value;
                    _sortedIndices[idx] = i;
                    idx++;
                }
            }

            // Sort by iRT (co-sort the index array)
            Array.Sort(_sortedIrts, _sortedIndices);
        }

        /// <summary>
        /// Returns the range of original precursor indices whose iRT falls within [irtLower, irtUpper].
        /// Uses binary search for O(log N) bounds finding.
        /// 
        /// The returned span of indices points into the original precursor list passed to the constructor.
        /// </summary>
        /// <param name="irtLower">Lower bound of the iRT window (inclusive).</param>
        /// <param name="irtUpper">Upper bound of the iRT window (inclusive).</param>
        /// <param name="startIndex">Output: first position in the sorted array within bounds.</param>
        /// <param name="count">Output: number of entries within bounds.</param>
        public void QueryRange(double irtLower, double irtUpper, out int startIndex, out int count)
        {
            startIndex = LowerBound(_sortedIrts, irtLower);
            int endIndex = UpperBound(_sortedIrts, irtUpper);
            count = endIndex - startIndex;
            if (count < 0) count = 0;
        }

        /// <summary>
        /// Returns original precursor indices for all entries within the iRT range.
        /// Allocates a new array; for hot paths, prefer QueryRange + GetOriginalIndex.
        /// </summary>
        public int[] GetOriginalIndices(double irtLower, double irtUpper)
        {
            QueryRange(irtLower, irtUpper, out int start, out int count);
            if (count == 0) return Array.Empty<int>();

            var result = new int[count];
            Array.Copy(_sortedIndices, start, result, 0, count);
            return result;
        }

        /// <summary>
        /// Returns the original precursor index for the i-th entry in sorted order.
        /// Valid for i in [startIndex, startIndex + count) from QueryRange.
        /// </summary>
        public int GetOriginalIndex(int sortedPosition) => _sortedIndices[sortedPosition];

        /// <summary>
        /// Returns the iRT value for the i-th entry in sorted order.
        /// </summary>
        public double GetIrt(int sortedPosition) => _sortedIrts[sortedPosition];

        /// <summary>
        /// Provides read-only access to the sorted iRT array (for diagnostics/plotting).
        /// </summary>
        public ReadOnlySpan<double> SortedIrts => _sortedIrts.AsSpan();

        // ── Binary search helpers ───────────────────────────────────────────

        /// <summary>
        /// Returns the index of the first element >= value (lower bound).
        /// If all elements are less than value, returns array.Length.
        /// </summary>
        private static int LowerBound(double[] sortedArray, double value)
        {
            int lo = 0, hi = sortedArray.Length;
            while (lo < hi)
            {
                int mid = lo + (hi - lo) / 2;
                if (sortedArray[mid] < value)
                    lo = mid + 1;
                else
                    hi = mid;
            }
            return lo;
        }

        /// <summary>
        /// Returns the index of the first element > value (upper bound).
        /// This is the exclusive end of the range [LowerBound, UpperBound).
        /// </summary>
        private static int UpperBound(double[] sortedArray, double value)
        {
            int lo = 0, hi = sortedArray.Length;
            while (lo < hi)
            {
                int mid = lo + (hi - lo) / 2;
                if (sortedArray[mid] <= value)
                    lo = mid + 1;
                else
                    hi = mid;
            }
            return lo;
        }
    }
}
