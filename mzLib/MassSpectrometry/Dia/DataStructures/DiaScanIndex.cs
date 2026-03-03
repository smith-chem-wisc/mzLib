// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Buffers;
using System.Collections.Generic;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Structure-of-Arrays (SoA) container for DIA MS2 scan data.
    /// 
    /// Instead of storing each scan as a separate object with its own m/z and intensity arrays
    /// (which causes poor cache locality and high GC pressure), this container packs all peak
    /// data into two contiguous float arrays (AllMz, AllIntensity). Each scan is then described
    /// by an offset and length into these arrays, plus its retention time and window ID.
    /// 
    /// This layout enables:
    ///   - Sequential memory access patterns for scan iteration
    ///   - Efficient binary search within a scan via Span slicing (no allocation)
    ///   - Direct GPU memory transfer (contiguous float arrays)
    ///   - Minimal GC pressure (a handful of arrays, not millions of objects)
    /// 
    /// Usage pattern:
    ///   1. Build via DiaScanIndexBuilder (single-pass from MsDataFile)
    ///   2. Access scan data via GetScanMzSpan / GetScanIntensitySpan
    ///   3. Get all scans for a window via GetScanIndicesForWindow
    /// </summary>
    public sealed class DiaScanIndex : IDisposable
    {
        // ── Contiguous peak data ────────────────────────────────────────────
        // All m/z values from all MS2 scans, packed end-to-end.
        // Sorted within each scan (as they come from the instrument).
        private readonly float[] _allMz;

        // All intensity values, aligned 1:1 with _allMz.
        private readonly float[] _allIntensity;

        // ── Per-scan metadata (parallel arrays, indexed by zero-based scan index) ──
        // Offset into _allMz/_allIntensity where this scan's peaks begin.
        private readonly int[] _scanOffsets;

        // Number of peaks in this scan.
        private readonly int[] _scanLengths;

        // Integer window ID for this scan (assigned during ingestion).
        private readonly int[] _scanWindowIds;

        // Retention time in minutes.
        private readonly float[] _scanRts;

        // One-based scan number from the original file (for traceability).
        private readonly int[] _scanOneBasedScanNumbers;

        // ── Window-to-scan index ────────────────────────────────────────────
        // Maps windowId → range of scan indices (start, count) in the per-scan arrays.
        // Scans within each window are sorted by RT.
        private readonly Dictionary<int, (int Start, int Count)> _windowToScanRange;

        // ── Isolation window definitions ────────────────────────────────────
        // Parallel arrays indexed by window ID.
        private readonly float[] _windowLowerBounds;
        private readonly float[] _windowUpperBounds;

        private bool _disposed;

        /// <summary>Total number of peaks across all scans.</summary>
        public int TotalPeakCount => _allMz.Length;

        /// <summary>Total number of MS2 scans in this index.</summary>
        public int ScanCount => _scanOffsets.Length;

        /// <summary>Number of distinct isolation windows.</summary>
        public int WindowCount => _windowLowerBounds.Length;

        /// <summary>
        /// Creates a DiaScanIndex from pre-built arrays. 
        /// Called by DiaScanIndexBuilder; not intended for direct use.
        /// All arrays are owned by this instance after construction.
        /// </summary>
        internal DiaScanIndex(
            float[] allMz,
            float[] allIntensity,
            int[] scanOffsets,
            int[] scanLengths,
            int[] scanWindowIds,
            float[] scanRts,
            int[] scanOneBasedScanNumbers,
            Dictionary<int, (int Start, int Count)> windowToScanRange,
            float[] windowLowerBounds,
            float[] windowUpperBounds)
        {
            _allMz = allMz ?? throw new ArgumentNullException(nameof(allMz));
            _allIntensity = allIntensity ?? throw new ArgumentNullException(nameof(allIntensity));
            _scanOffsets = scanOffsets ?? throw new ArgumentNullException(nameof(scanOffsets));
            _scanLengths = scanLengths ?? throw new ArgumentNullException(nameof(scanLengths));
            _scanWindowIds = scanWindowIds ?? throw new ArgumentNullException(nameof(scanWindowIds));
            _scanRts = scanRts ?? throw new ArgumentNullException(nameof(scanRts));
            _scanOneBasedScanNumbers = scanOneBasedScanNumbers ?? throw new ArgumentNullException(nameof(scanOneBasedScanNumbers));
            _windowToScanRange = windowToScanRange ?? throw new ArgumentNullException(nameof(windowToScanRange));
            _windowLowerBounds = windowLowerBounds ?? throw new ArgumentNullException(nameof(windowLowerBounds));
            _windowUpperBounds = windowUpperBounds ?? throw new ArgumentNullException(nameof(windowUpperBounds));

            if (allMz.Length != allIntensity.Length)
                throw new ArgumentException("m/z and intensity arrays must have the same length.");
            if (scanOffsets.Length != scanLengths.Length ||
                scanOffsets.Length != scanWindowIds.Length ||
                scanOffsets.Length != scanRts.Length ||
                scanOffsets.Length != scanOneBasedScanNumbers.Length)
                throw new ArgumentException("All per-scan arrays must have the same length.");
            if (windowLowerBounds.Length != windowUpperBounds.Length)
                throw new ArgumentException("Window bound arrays must have the same length.");
        }

        // ── Scan-level accessors ────────────────────────────────────────────

        /// <summary>
        /// Returns the m/z values for the specified scan as a read-only span.
        /// Zero-copy: this is a view into the contiguous array, not a new allocation.
        /// </summary>
        public ReadOnlySpan<float> GetScanMzSpan(int scanIndex)
        {
            return new ReadOnlySpan<float>(_allMz, _scanOffsets[scanIndex], _scanLengths[scanIndex]);
        }

        /// <summary>
        /// Returns the intensity values for the specified scan as a read-only span.
        /// Zero-copy: this is a view into the contiguous array, not a new allocation.
        /// </summary>
        public ReadOnlySpan<float> GetScanIntensitySpan(int scanIndex)
        {
            return new ReadOnlySpan<float>(_allIntensity, _scanOffsets[scanIndex], _scanLengths[scanIndex]);
        }

        /// <summary>Returns the retention time (minutes) for the specified scan.</summary>
        public float GetScanRt(int scanIndex) => _scanRts[scanIndex];

        /// <summary>Returns the window ID for the specified scan.</summary>
        public int GetScanWindowId(int scanIndex) => _scanWindowIds[scanIndex];

        /// <summary>Returns the original one-based scan number for traceability.</summary>
        public int GetScanOneBasedScanNumber(int scanIndex) => _scanOneBasedScanNumbers[scanIndex];

        /// <summary>Returns the number of peaks in the specified scan.</summary>
        public int GetScanPeakCount(int scanIndex) => _scanLengths[scanIndex];

        // ── Window-level accessors ──────────────────────────────────────────

        /// <summary>
        /// Returns the range of zero-based scan indices that belong to the given window ID.
        /// Scans within this range are sorted by retention time.
        /// Returns false if the window ID is not found.
        /// </summary>
        public bool TryGetScanRangeForWindow(int windowId, out int startScanIndex, out int scanCount)
        {
            if (_windowToScanRange.TryGetValue(windowId, out var range))
            {
                startScanIndex = range.Start;
                scanCount = range.Count;
                return true;
            }
            startScanIndex = 0;
            scanCount = 0;
            return false;
        }

        /// <summary>
        /// Returns the isolation window m/z boundaries for the given window ID.
        /// </summary>
        public (float LowerBound, float UpperBound) GetWindowBounds(int windowId)
        {
            if (windowId < 0 || windowId >= _windowLowerBounds.Length)
                throw new ArgumentOutOfRangeException(nameof(windowId));
            return (_windowLowerBounds[windowId], _windowUpperBounds[windowId]);
        }

        /// <summary>Returns all window IDs present in this index.</summary>
        public IReadOnlyCollection<int> GetWindowIds() => _windowToScanRange.Keys;

        // ── Bulk data access (for GPU transfer or batch operations) ─────────

        /// <summary>
        /// Provides direct read-only access to the entire contiguous m/z array.
        /// Intended for GPU memory transfer or SIMD batch operations.
        /// </summary>
        public ReadOnlyMemory<float> AllMz => _allMz.AsMemory();

        /// <summary>
        /// Provides direct read-only access to the entire contiguous intensity array.
        /// Intended for GPU memory transfer or SIMD batch operations.
        /// </summary>
        public ReadOnlyMemory<float> AllIntensity => _allIntensity.AsMemory();

        /// <summary>
        /// Provides direct read-only access to scan retention times.
        /// </summary>
        public ReadOnlySpan<float> AllScanRts => _scanRts.AsSpan();

        /// <summary>
        /// Provides direct read-only access to scan window IDs.
        /// </summary>
        public ReadOnlySpan<int> AllScanWindowIds => _scanWindowIds.AsSpan();

        /// <summary>
        /// Provides direct read-only access to per-scan peak offsets into AllMz/AllIntensity.
        /// Intended for GPU memory transfer where the full offset array must be uploaded.
        /// </summary>
        public ReadOnlySpan<int> AllScanOffsets => _scanOffsets.AsSpan();

        /// <summary>
        /// Provides direct read-only access to per-scan peak counts.
        /// Intended for GPU memory transfer.
        /// </summary>
        public ReadOnlySpan<int> AllScanLengths => _scanLengths.AsSpan();

        /// <summary>
        /// Finds the DIA isolation window that contains the given precursor m/z.
        /// Linear scan over window bounds; O(W) where W is typically 20-60.
        /// Returns -1 if the precursor falls outside all windows.
        /// </summary>
        public int FindWindowForPrecursorMz(float precursorMz)
        {
            for (int w = 0; w < _windowLowerBounds.Length; w++)
            {
                if (precursorMz >= _windowLowerBounds[w] && precursorMz <= _windowUpperBounds[w])
                    return w;
            }
            return -1;
        }

        /// <summary>
        /// Overload accepting double for convenience (e.g., LibrarySpectrum.PrecursorMz).
        /// </summary>
        public int FindWindowForPrecursorMz(double precursorMz)
        {
            return FindWindowForPrecursorMz((float)precursorMz);
        }

        /// <summary>
        /// Returns the global minimum RT across all scans, or 0 if empty.
        /// Scans are sorted by (windowId, RT), so the first scan's RT is the global min.
        /// </summary>
        public float GetGlobalRtMin()
        {
            if (ScanCount == 0) return 0f;
            float min = _scanRts[0];
            for (int i = 1; i < _scanRts.Length; i++)
                if (_scanRts[i] < min) min = _scanRts[i];
            return min;
        }

        /// <summary>
        /// Returns the global maximum RT across all scans, or 0 if empty.
        /// </summary>
        public float GetGlobalRtMax()
        {
            if (ScanCount == 0) return 0f;
            float max = _scanRts[0];
            for (int i = 1; i < _scanRts.Length; i++)
                if (_scanRts[i] > max) max = _scanRts[i];
            return max;
        }

        public void Dispose()
        {
            // Currently no pooled arrays to return.
            // This is here for forward-compatibility when we switch to ArrayPool<T>.
            _disposed = true;
        }
    }
}
