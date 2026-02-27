// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Runtime.CompilerServices;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// CPU implementation of fragment ion extraction from DIA data.
    /// 
    /// For each FragmentQuery, this extractor:
    ///   1. Looks up the scan range for the query's window ID
    ///   2. Iterates only scans within the query's RT window
    ///   3. For each qualifying scan, uses binary search to find the m/z region 
    ///      within TargetMz ± ppm tolerance
    ///   4. Sums intensities of all matching peaks and writes the best (highest 
    ///      intensity) match as the XIC data point for that scan
    /// 
    /// Performance characteristics:
    ///   - O(S × log P) per query, where S = scans in RT window, P = peaks per scan
    ///   - No heap allocations in the extraction loop
    ///   - All data access is via Span slicing into the contiguous SoA arrays
    ///   - Suitable for Parallel.For at the window level (one instance per thread)
    /// 
    /// Thread safety: NOT thread-safe. Use one instance per thread, or synchronize externally.
    /// The DiaScanIndex it reads from is immutable and safe to share across threads.
    /// </summary>
    public sealed class CpuFragmentExtractor : IFragmentExtractor
    {
        private readonly DiaScanIndex _index;
        private bool _disposed;

        /// <summary>
        /// Creates a CPU fragment extractor backed by the given DIA scan index.
        /// The index must remain alive for the lifetime of this extractor.
        /// </summary>
        public CpuFragmentExtractor(DiaScanIndex index)
        {
            _index = index ?? throw new ArgumentNullException(nameof(index));
        }

        /// <summary>
        /// Extracts fragment ion chromatograms for a batch of queries.
        /// 
        /// For each query:
        ///   - Finds the scan range for the query's window ID
        ///   - Skips scans outside [RtMin, RtMax]
        ///   - Binary searches each qualifying scan's m/z array for peaks within tolerance
        ///   - For each scan with a match, writes one XIC data point (RT, summed intensity)
        ///     into the output buffers
        /// 
        /// Returns the total number of XIC data points written across all queries.
        /// The caller must ensure rtBuffer and intensityBuffer are large enough.
        /// A safe size is: queries.Length × maxScansPerWindow.
        /// </summary>
        public int ExtractBatch(
            ReadOnlySpan<FragmentQuery> queries,
            Span<FragmentResult> results,
            Span<float> rtBuffer,
            Span<float> intensityBuffer)
        {
            if (queries.Length != results.Length)
                throw new ArgumentException("Queries and results spans must have the same length.");

            int totalDataPoints = 0;

            for (int q = 0; q < queries.Length; q++)
            {
                ref readonly var query = ref queries[q];
                int dataPointStart = totalDataPoints;
                float totalIntensity = 0f;

                // Look up scan range for this window
                if (!_index.TryGetScanRangeForWindow(query.WindowId, out int scanStart, out int scanCount))
                {
                    // Unknown window — no data
                    results[q] = new FragmentResult(query.QueryId, 0, dataPointStart, dataPointStart, 0f);
                    continue;
                }

                // Compute m/z bounds from ppm tolerance
                float mzLow = PpmToMzLow(query.TargetMz, query.TolerancePpm);
                float mzHigh = PpmToMzHigh(query.TargetMz, query.TolerancePpm);

                int scanEnd = scanStart + scanCount;
                for (int s = scanStart; s < scanEnd; s++)
                {
                    // RT filtering: scans are sorted by RT within each window,
                    // so we can skip early and break early.
                    float rt = _index.GetScanRt(s);
                    if (rt < query.RtMin) continue;
                    if (rt > query.RtMax) break;

                    // Binary search for matching peaks in this scan
                    ReadOnlySpan<float> mzSpan = _index.GetScanMzSpan(s);
                    if (mzSpan.Length == 0) continue;

                    // Find the first index where mz >= mzLow
                    int lo = LowerBound(mzSpan, mzLow);
                    if (lo >= mzSpan.Length) continue;
                    if (mzSpan[lo] > mzHigh) continue;

                    // Sum intensities of all peaks in [mzLow, mzHigh]
                    ReadOnlySpan<float> intSpan = _index.GetScanIntensitySpan(s);
                    float scanIntensity = 0f;
                    for (int p = lo; p < mzSpan.Length && mzSpan[p] <= mzHigh; p++)
                    {
                        scanIntensity += intSpan[p];
                    }

                    if (scanIntensity > 0f)
                    {
                        rtBuffer[totalDataPoints] = rt;
                        intensityBuffer[totalDataPoints] = scanIntensity;
                        totalDataPoints++;
                        totalIntensity += scanIntensity;
                    }
                }

                int dataPointCount = totalDataPoints - dataPointStart;
                results[q] = new FragmentResult(
                    query.QueryId, dataPointCount, dataPointStart, dataPointStart, totalIntensity);
            }

            return totalDataPoints;
        }

        /// <summary>
        /// Binary search for the first index in a sorted span where value >= target.
        /// Returns span.Length if all values are less than target.
        /// 
        /// This is a standard lower-bound binary search with no allocations.
        /// Used to find the start of the m/z range that could match a fragment query.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static int LowerBound(ReadOnlySpan<float> sorted, float target)
        {
            int lo = 0;
            int hi = sorted.Length;
            while (lo < hi)
            {
                int mid = lo + ((hi - lo) >> 1);
                if (sorted[mid] < target)
                    lo = mid + 1;
                else
                    hi = mid;
            }
            return lo;
        }

        /// <summary>
        /// Converts a target m/z and ppm tolerance to the lower m/z bound.
        /// Formula: mz × (1 - ppm / 1e6)
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static float PpmToMzLow(float mz, float ppm)
        {
            return mz * (1f - ppm / 1_000_000f);
        }

        /// <summary>
        /// Converts a target m/z and ppm tolerance to the upper m/z bound.
        /// Formula: mz × (1 + ppm / 1e6)
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static float PpmToMzHigh(float mz, float ppm)
        {
            return mz * (1f + ppm / 1_000_000f);
        }

        public void Dispose()
        {
            _disposed = true;
            // No unmanaged resources. DiaScanIndex is not owned by this extractor.
        }
    }
}
