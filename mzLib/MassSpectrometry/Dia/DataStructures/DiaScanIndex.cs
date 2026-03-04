// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Buffers;
using System.Collections.Generic;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Structure-of-Arrays (SoA) container for DIA MS2 scan data AND MS1 scan data.
    /// 
    /// Instead of storing each scan as a separate object with its own m/z and intensity arrays
    /// (which causes poor cache locality and high GC pressure), this container packs all peak
    /// data into contiguous float arrays. Each scan is then described by an offset and length
    /// into these arrays, plus its retention time.
    /// 
    /// MS2 layout:
    ///   - AllMz / AllIntensity: contiguous peak data for all MS2 scans
    ///   - Per-scan: offset, length, window ID, RT, one-based scan number
    ///   - Window-to-scan dictionary for O(1) window lookup
    ///   - Scans within each window are sorted by RT
    /// 
    /// MS1 layout (Phase 16B, Prompt 3):
    ///   - _ms1AllMz / _ms1AllIntensity: contiguous peak data for all MS1 scans
    ///   - Per-scan: offset, length, RT, one-based scan number
    ///   - MS1 scans are sorted globally by RT (no window association)
    ///   - FindMs1ScanIndexAtRt performs binary search for O(log N) RT lookup
    /// 
    /// This layout enables:
    ///   - Sequential memory access patterns for scan iteration
    ///   - Efficient binary search within a scan via Span slicing (no allocation)
    ///   - Direct GPU memory transfer (contiguous float arrays)
    ///   - Minimal GC pressure (a handful of arrays, not millions of objects)
    /// 
    /// Usage pattern:
    ///   1. Build via DiaScanIndexBuilder (single-pass from MsDataFile)
    ///   2. Access MS2 scan data via GetScanMzSpan / GetScanIntensitySpan
    ///   3. Access MS1 scan data via GetMs1ScanPeaks / FindMs1ScanIndexAtRt
    ///   4. Get all MS2 scans for a window via GetScanIndicesForWindow
    /// </summary>
    public sealed class DiaScanIndex : IDisposable
    {
        // ── MS2: Contiguous peak data ───────────────────────────────────────
        // All m/z values from all MS2 scans, packed end-to-end.
        // Sorted within each scan (as they come from the instrument).
        private readonly float[] _allMz;

        // All intensity values, aligned 1:1 with _allMz.
        private readonly float[] _allIntensity;

        // ── MS2: Per-scan metadata (parallel arrays, indexed by zero-based scan index) ──
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

        // ── MS2: Window-to-scan index ───────────────────────────────────────
        // Maps windowId → range of scan indices (start, count) in the per-scan arrays.
        // Scans within each window are sorted by RT.
        private readonly Dictionary<int, (int Start, int Count)> _windowToScanRange;

        // ── MS2: Isolation window definitions ──────────────────────────────
        // Parallel arrays indexed by window ID.
        private readonly float[] _windowLowerBounds;
        private readonly float[] _windowUpperBounds;

        // ── MS1: Contiguous peak data ───────────────────────────────────────
        // All m/z values from all MS1 scans, packed end-to-end.
        // Sorted within each scan (as they come from the instrument).
        private readonly float[] _ms1AllMz;

        // All intensity values, aligned 1:1 with _ms1AllMz.
        private readonly float[] _ms1AllIntensity;

        // ── MS1: Per-scan metadata (parallel arrays, indexed by zero-based MS1 scan index) ──
        // Offset into _ms1AllMz/_ms1AllIntensity where this scan's peaks begin.
        private readonly int[] _ms1ScanOffsets;

        // Number of peaks in this MS1 scan.
        private readonly int[] _ms1ScanLengths;

        // Retention time in minutes. MS1 scans are sorted ascending by RT,
        // enabling binary search in FindMs1ScanIndexAtRt.
        private readonly float[] _ms1ScanRts;

        // One-based scan number from the original file (for traceability).
        private readonly int[] _ms1OneBasedScanNumbers;

        private bool _disposed;

        // ── MS2: Public properties ──────────────────────────────────────────

        /// <summary>Total number of peaks across all MS2 scans.</summary>
        public int TotalPeakCount => _allMz.Length;

        /// <summary>Total number of MS2 scans in this index.</summary>
        public int ScanCount => _scanOffsets.Length;

        /// <summary>Number of distinct isolation windows.</summary>
        public int WindowCount => _windowLowerBounds.Length;

        // ── MS1: Public properties ──────────────────────────────────────────

        /// <summary>Total number of MS1 scans in this index.</summary>
        public int Ms1ScanCount => _ms1ScanOffsets.Length;

        /// <summary>Total number of peaks across all MS1 scans.</summary>
        public int Ms1TotalPeakCount => _ms1AllMz.Length;

        /// <summary>
        /// Creates a DiaScanIndex from pre-built arrays.
        /// Called by DiaScanIndexBuilder; not intended for direct use.
        /// All arrays are owned by this instance after construction.
        /// </summary>
        internal DiaScanIndex(
            // MS2 arrays
            float[] allMz,
            float[] allIntensity,
            int[] scanOffsets,
            int[] scanLengths,
            int[] scanWindowIds,
            float[] scanRts,
            int[] scanOneBasedScanNumbers,
            Dictionary<int, (int Start, int Count)> windowToScanRange,
            float[] windowLowerBounds,
            float[] windowUpperBounds,
            // MS1 arrays
            float[] ms1AllMz,
            float[] ms1AllIntensity,
            int[] ms1ScanOffsets,
            int[] ms1ScanLengths,
            float[] ms1ScanRts,
            int[] ms1OneBasedScanNumbers)
        {
            // MS2 validation
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
                throw new ArgumentException("MS2 m/z and intensity arrays must have the same length.");
            if (scanOffsets.Length != scanLengths.Length ||
                scanOffsets.Length != scanWindowIds.Length ||
                scanOffsets.Length != scanRts.Length ||
                scanOffsets.Length != scanOneBasedScanNumbers.Length)
                throw new ArgumentException("All MS2 per-scan arrays must have the same length.");
            if (windowLowerBounds.Length != windowUpperBounds.Length)
                throw new ArgumentException("Window bound arrays must have the same length.");

            // MS1 validation
            _ms1AllMz = ms1AllMz ?? throw new ArgumentNullException(nameof(ms1AllMz));
            _ms1AllIntensity = ms1AllIntensity ?? throw new ArgumentNullException(nameof(ms1AllIntensity));
            _ms1ScanOffsets = ms1ScanOffsets ?? throw new ArgumentNullException(nameof(ms1ScanOffsets));
            _ms1ScanLengths = ms1ScanLengths ?? throw new ArgumentNullException(nameof(ms1ScanLengths));
            _ms1ScanRts = ms1ScanRts ?? throw new ArgumentNullException(nameof(ms1ScanRts));
            _ms1OneBasedScanNumbers = ms1OneBasedScanNumbers ?? throw new ArgumentNullException(nameof(ms1OneBasedScanNumbers));

            if (ms1AllMz.Length != ms1AllIntensity.Length)
                throw new ArgumentException("MS1 m/z and intensity arrays must have the same length.");
            if (ms1ScanOffsets.Length != ms1ScanLengths.Length ||
                ms1ScanOffsets.Length != ms1ScanRts.Length ||
                ms1ScanOffsets.Length != ms1OneBasedScanNumbers.Length)
                throw new ArgumentException("All MS1 per-scan arrays must have the same length.");

            // Validate MS1 scans are sorted by RT (required for binary search)
            for (int i = 1; i < ms1ScanRts.Length; i++)
            {
                if (ms1ScanRts[i] < ms1ScanRts[i - 1])
                    throw new ArgumentException(
                        $"MS1 scans must be sorted ascending by RT. " +
                        $"Scan at index {i} has RT={ms1ScanRts[i]:F4} < previous RT={ms1ScanRts[i - 1]:F4}. " +
                        "Sort MS1 scans by RetentionTime before calling DiaScanIndexBuilder.Build().");
            }
        }

        // ── MS2: Scan-level accessors ───────────────────────────────────────

        /// <summary>
        /// Returns the m/z values for the specified MS2 scan as a read-only span.
        /// Zero-copy: this is a view into the contiguous array, not a new allocation.
        /// </summary>
        public ReadOnlySpan<float> GetScanMzSpan(int scanIndex)
        {
            return new ReadOnlySpan<float>(_allMz, _scanOffsets[scanIndex], _scanLengths[scanIndex]);
        }

        /// <summary>
        /// Returns the intensity values for the specified MS2 scan as a read-only span.
        /// Zero-copy: this is a view into the contiguous array, not a new allocation.
        /// </summary>
        public ReadOnlySpan<float> GetScanIntensitySpan(int scanIndex)
        {
            return new ReadOnlySpan<float>(_allIntensity, _scanOffsets[scanIndex], _scanLengths[scanIndex]);
        }

        /// <summary>Returns the retention time (minutes) for the specified MS2 scan.</summary>
        public float GetScanRt(int scanIndex) => _scanRts[scanIndex];

        /// <summary>Returns the window ID for the specified MS2 scan.</summary>
        public int GetScanWindowId(int scanIndex) => _scanWindowIds[scanIndex];

        /// <summary>Returns the original one-based scan number for the specified MS2 scan.</summary>
        public int GetScanOneBasedScanNumber(int scanIndex) => _scanOneBasedScanNumbers[scanIndex];

        /// <summary>Returns the number of peaks in the specified MS2 scan.</summary>
        public int GetScanPeakCount(int scanIndex) => _scanLengths[scanIndex];

        // ── MS1: Scan-level accessors ───────────────────────────────────────

        /// <summary>
        /// Retrieves the m/z and intensity peak arrays for the specified MS1 scan.
        /// Zero-copy: both spans are views into contiguous arrays, not new allocations.
        /// </summary>
        /// <param name="ms1ScanIndex">Zero-based MS1 scan index.</param>
        /// <param name="mzs">Receives the m/z values for this scan.</param>
        /// <param name="intensities">Receives the intensity values for this scan.</param>
        public void GetMs1ScanPeaks(
            int ms1ScanIndex,
            out ReadOnlySpan<float> mzs,
            out ReadOnlySpan<float> intensities)
        {
            if ((uint)ms1ScanIndex >= (uint)_ms1ScanOffsets.Length)
                throw new ArgumentOutOfRangeException(nameof(ms1ScanIndex),
                    $"MS1 scan index {ms1ScanIndex} is out of range [0, {_ms1ScanOffsets.Length}).");

            int offset = _ms1ScanOffsets[ms1ScanIndex];
            int length = _ms1ScanLengths[ms1ScanIndex];
            mzs = new ReadOnlySpan<float>(_ms1AllMz, offset, length);
            intensities = new ReadOnlySpan<float>(_ms1AllIntensity, offset, length);
        }

        /// <summary>Returns the retention time (minutes) for the specified MS1 scan.</summary>
        /// <param name="ms1ScanIndex">Zero-based MS1 scan index.</param>
        public float GetMs1ScanRt(int ms1ScanIndex)
        {
            if ((uint)ms1ScanIndex >= (uint)_ms1ScanRts.Length)
                throw new ArgumentOutOfRangeException(nameof(ms1ScanIndex),
                    $"MS1 scan index {ms1ScanIndex} is out of range [0, {_ms1ScanRts.Length}).");
            return _ms1ScanRts[ms1ScanIndex];
        }

        /// <summary>Returns the number of peaks in the specified MS1 scan.</summary>
        public int GetMs1ScanPeakCount(int ms1ScanIndex) => _ms1ScanLengths[ms1ScanIndex];

        /// <summary>Returns the original one-based scan number for the specified MS1 scan.</summary>
        public int GetMs1ScanOneBasedScanNumber(int ms1ScanIndex) => _ms1OneBasedScanNumbers[ms1ScanIndex];

        /// <summary>
        /// Binary search: finds the index of the first MS1 scan with RT ≥ rtMin.
        /// 
        /// MS1 scans are stored sorted ascending by RT, so a lower_bound binary search
        /// directly locates the start of the RT window in O(log N).
        /// 
        /// Typical usage (iterate MS1 scans in [rtMin, rtMax]):
        /// <code>
        ///   int start = index.FindMs1ScanIndexAtRt(rtMin);
        ///   for (int i = start; i &lt; index.Ms1ScanCount; i++)
        ///   {
        ///       if (index.GetMs1ScanRt(i) > rtMax) break;
        ///       index.GetMs1ScanPeaks(i, out var mzs, out var intensities);
        ///       // ... process peaks ...
        ///   }
        /// </code>
        /// </summary>
        /// <param name="rtMin">Lower RT bound (minutes, inclusive).</param>
        /// <returns>
        /// Zero-based index of the first MS1 scan with RT ≥ rtMin.
        /// Returns Ms1ScanCount if all scans have RT &lt; rtMin (i.e., rtMin is past the run end).
        /// Returns 0 if rtMin ≤ the RT of the first scan.
        /// </returns>
        public int FindMs1ScanIndexAtRt(float rtMin)
        {
            if (_ms1ScanRts.Length == 0) return 0;

            // Standard lower_bound binary search
            int lo = 0;
            int hi = _ms1ScanRts.Length; // exclusive upper bound

            while (lo < hi)
            {
                int mid = lo + ((hi - lo) >> 1);
                if (_ms1ScanRts[mid] < rtMin)
                    lo = mid + 1;
                else
                    hi = mid;
            }

            return lo; // index of first element with RT >= rtMin
        }

        // ── MS2: Window-level accessors ─────────────────────────────────────

        /// <summary>
        /// Returns the range of zero-based MS2 scan indices that belong to the given window ID.
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

        // ── MS2: Bulk data access (for GPU transfer or batch operations) ────

        /// <summary>
        /// Provides direct read-only access to the entire contiguous MS2 m/z array.
        /// Intended for GPU memory transfer or SIMD batch operations.
        /// </summary>
        public ReadOnlyMemory<float> AllMz => _allMz.AsMemory();

        /// <summary>
        /// Provides direct read-only access to the entire contiguous MS2 intensity array.
        /// Intended for GPU memory transfer or SIMD batch operations.
        /// </summary>
        public ReadOnlyMemory<float> AllIntensity => _allIntensity.AsMemory();

        /// <summary>Provides direct read-only access to MS2 scan retention times.</summary>
        public ReadOnlySpan<float> AllScanRts => _scanRts.AsSpan();

        /// <summary>Provides direct read-only access to MS2 scan window IDs.</summary>
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

        // ── MS1: Bulk data access ───────────────────────────────────────────

        /// <summary>
        /// Provides direct read-only access to the entire contiguous MS1 m/z array.
        /// Intended for GPU memory transfer or SIMD batch operations.
        /// </summary>
        public ReadOnlyMemory<float> AllMs1Mz => _ms1AllMz.AsMemory();

        /// <summary>
        /// Provides direct read-only access to the entire contiguous MS1 intensity array.
        /// Intended for GPU memory transfer or SIMD batch operations.
        /// </summary>
        public ReadOnlyMemory<float> AllMs1Intensity => _ms1AllIntensity.AsMemory();

        /// <summary>
        /// Provides direct read-only access to MS1 scan retention times (sorted ascending).
        /// </summary>
        public ReadOnlySpan<float> AllMs1ScanRts => _ms1ScanRts.AsSpan();

        /// <summary>Provides direct read-only access to MS1 per-scan peak offsets.</summary>
        public ReadOnlySpan<int> AllMs1ScanOffsets => _ms1ScanOffsets.AsSpan();

        /// <summary>Provides direct read-only access to MS1 per-scan peak counts.</summary>
        public ReadOnlySpan<int> AllMs1ScanLengths => _ms1ScanLengths.AsSpan();

        // ── Utility methods ─────────────────────────────────────────────────

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
        /// Returns the global minimum RT across all MS2 scans, or 0 if empty.
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
        /// Returns the global maximum RT across all MS2 scans, or 0 if empty.
        /// </summary>
        public float GetGlobalRtMax()
        {
            if (ScanCount == 0) return 0f;
            float max = _scanRts[0];
            for (int i = 1; i < _scanRts.Length; i++)
                if (_scanRts[i] > max) max = _scanRts[i];
            return max;
        }

        /// <summary>
        /// Returns the global minimum RT across all MS1 scans, or 0 if empty.
        /// Since MS1 scans are sorted by RT, this is simply the first scan's RT.
        /// </summary>
        public float GetMs1GlobalRtMin()
        {
            return Ms1ScanCount == 0 ? 0f : _ms1ScanRts[0];
        }

        /// <summary>
        /// Returns the global maximum RT across all MS1 scans, or 0 if empty.
        /// Since MS1 scans are sorted by RT, this is simply the last scan's RT.
        /// </summary>
        public float GetMs1GlobalRtMax()
        {
            return Ms1ScanCount == 0 ? 0f : _ms1ScanRts[Ms1ScanCount - 1];
        }

        public void Dispose()
        {
            // Currently no pooled arrays to return.
            // This is here for forward-compatibility when we switch to ArrayPool<T>.
            _disposed = true;
        }
    }
}