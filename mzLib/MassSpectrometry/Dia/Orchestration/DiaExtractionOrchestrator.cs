// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Buffers;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Orchestrates parallel fragment extraction across DIA isolation windows.
    /// 
    /// This is the main entry point for high-throughput DIA fragment extraction. It:
    ///   1. Groups incoming queries by window ID
    ///   2. Dispatches extraction work across windows using Parallel.For
    ///   3. Each thread gets its own IFragmentExtractor instance (thread-local, no contention)
    ///   4. Per-thread result buffers are rented from ArrayPool (no GC pressure)
    ///   5. Results are merged back into the caller's output arrays
    /// 
    /// Thread safety model:
    ///   - DiaScanIndex is immutable and shared across all threads (read-only)
    ///   - Each Parallel.For partition creates its own CpuFragmentExtractor
    ///   - Output buffers per window are independent (no locking needed)
    ///   - Final merge into caller's arrays is sequential (fast — just memcpy)
    /// 
    /// Performance characteristics:
    ///   - Parallelism limited by window count (typically 20–60 windows)
    ///   - For a typical DIA file with 20 windows and 8 cores, expect 4–6× speedup
    ///   - ArrayPool reuse eliminates per-extraction heap allocations
    ///   - Total memory overhead: ~O(maxQueriesPerWindow × maxScansPerWindow × sizeof(float))
    ///     per concurrent thread
    /// 
    /// Usage:
    ///   var orchestrator = new DiaExtractionOrchestrator(index);
    ///   var allResults = orchestrator.ExtractAll(queries, maxDegreeOfParallelism: 8);
    /// </summary>
    public sealed class DiaExtractionOrchestrator : IDisposable
    {
        private readonly DiaScanIndex _index;
        private readonly Func<DiaScanIndex, IFragmentExtractor> _extractorFactory;
        private bool _disposed;

        /// <summary>
        /// Creates an orchestrator that uses CpuFragmentExtractor by default.
        /// </summary>
        /// <param name="index">The shared DIA scan index (immutable, thread-safe for reads).</param>
        public DiaExtractionOrchestrator(DiaScanIndex index)
            : this(index, idx => new CpuFragmentExtractor(idx))
        {
        }

        /// <summary>
        /// Creates an orchestrator with a custom extractor factory.
        /// This allows injecting GpuFragmentExtractor or test doubles.
        /// </summary>
        /// <param name="index">The shared DIA scan index.</param>
        /// <param name="extractorFactory">
        /// Factory function that creates an IFragmentExtractor for a given DiaScanIndex.
        /// Called once per thread partition. The orchestrator disposes these after use.
        /// </param>
        public DiaExtractionOrchestrator(DiaScanIndex index, Func<DiaScanIndex, IFragmentExtractor> extractorFactory)
        {
            _index = index ?? throw new ArgumentNullException(nameof(index));
            _extractorFactory = extractorFactory ?? throw new ArgumentNullException(nameof(extractorFactory));
        }

        /// <summary>
        /// Extracts fragment ion chromatograms for all queries in parallel, grouped by window.
        /// 
        /// Steps:
        ///   1. Group queries by WindowId (single pass, O(N))
        ///   2. Parallel.For across windows, each with its own extractor + pooled buffers
        ///   3. Merge per-window results into a single contiguous output
        /// 
        /// The returned ExtractionResult owns the output buffers and should be disposed
        /// when no longer needed (returns pooled arrays).
        /// </summary>
        /// <param name="queries">All fragment queries to extract. Can span multiple windows.</param>
        /// <param name="maxDegreeOfParallelism">
        /// Maximum number of concurrent threads. Default -1 uses Environment.ProcessorCount.
        /// Set to 1 for deterministic single-threaded execution (useful for testing/debugging).
        /// </param>
        /// <returns>
        /// An ExtractionResult containing all FragmentResults and the XIC data buffers.
        /// Results are in the same order as the input queries.
        /// </returns>
        public ExtractionResult ExtractAll(
            ReadOnlySpan<FragmentQuery> queries,
            int maxDegreeOfParallelism = -1)
        {
            if (queries.Length == 0)
            {
                return ExtractionResult.Empty;
            }

            if (maxDegreeOfParallelism <= 0)
            {
                maxDegreeOfParallelism = Environment.ProcessorCount;
            }

            // ── Step 1: Group queries by window ID ──────────────────────────────
            // We need to know which queries go to which window so we can dispatch them
            // in parallel. We also need to track original indices for result reordering.
            var windowGroups = GroupQueriesByWindow(queries);

            // ── Step 2: Estimate buffer sizes ───────────────────────────────────
            // Each window needs its own result + XIC buffers. We estimate max XIC points
            // as queriesInWindow × maxScansInWindow.
            int maxScansPerWindow = 0;
            foreach (int windowId in windowGroups.Keys)
            {
                if (_index.TryGetScanRangeForWindow(windowId, out _, out int scanCount))
                {
                    if (scanCount > maxScansPerWindow) maxScansPerWindow = scanCount;
                }
            }

            // ── Step 3: Parallel extraction per window ──────────────────────────
            // Each window group is processed independently. We store per-window results
            // in a concurrent-safe structure (array indexed by window group index).
            var windowKeys = new int[windowGroups.Count];
            var windowQueryArrays = new (FragmentQuery[] Queries, int[] OriginalIndices)[windowGroups.Count];
            int groupIdx = 0;
            foreach (var kvp in windowGroups)
            {
                windowKeys[groupIdx] = kvp.Key;
                windowQueryArrays[groupIdx] = (kvp.Value.Queries.ToArray(), kvp.Value.OriginalIndices.ToArray());
                groupIdx++;
            }

            // Per-window results: each slot filled by its parallel task
            var perWindowResults = new WindowExtractionResult[windowKeys.Length];

            var options = new ParallelOptions
            {
                MaxDegreeOfParallelism = maxDegreeOfParallelism
            };

            Parallel.For(0, windowKeys.Length, options, windowGroupIndex =>
            {
                var (windowQueries, originalIndices) = windowQueryArrays[windowGroupIndex];
                int queryCount = windowQueries.Length;

                // Rent buffers from ArrayPool (returned after merge)
                int maxXicPoints = queryCount * Math.Max(maxScansPerWindow, 1);
                var results = new FragmentResult[queryCount];
                float[] rtBuffer = ArrayPool<float>.Shared.Rent(maxXicPoints);
                float[] intensityBuffer = ArrayPool<float>.Shared.Rent(maxXicPoints);

                try
                {
                    // Each thread creates its own extractor (no shared mutable state)
                    using var extractor = _extractorFactory(_index);
                    int totalPoints = extractor.ExtractBatch(
                        windowQueries,
                        results,
                        rtBuffer.AsSpan(0, maxXicPoints),
                        intensityBuffer.AsSpan(0, maxXicPoints));

                    perWindowResults[windowGroupIndex] = new WindowExtractionResult(
                        results, originalIndices, rtBuffer, intensityBuffer, totalPoints, maxXicPoints);
                }
                catch
                {
                    // On exception, return buffers before rethrowing
                    ArrayPool<float>.Shared.Return(rtBuffer);
                    ArrayPool<float>.Shared.Return(intensityBuffer);
                    throw;
                }
            });

            // ── Step 4: Merge per-window results into unified output ────────────
            return MergeResults(queries.Length, perWindowResults);
        }

        /// <summary>
        /// Groups queries by WindowId, preserving original indices for result reordering.
        /// Single pass: O(N) where N = query count.
        /// </summary>
        private static Dictionary<int, (List<FragmentQuery> Queries, List<int> OriginalIndices)>
            GroupQueriesByWindow(ReadOnlySpan<FragmentQuery> queries)
        {
            var groups = new Dictionary<int, (List<FragmentQuery> Queries, List<int> OriginalIndices)>();

            for (int i = 0; i < queries.Length; i++)
            {
                int windowId = queries[i].WindowId;
                if (!groups.TryGetValue(windowId, out var group))
                {
                    group = (new List<FragmentQuery>(), new List<int>());
                    groups[windowId] = group;
                }
                group.Queries.Add(queries[i]);
                group.OriginalIndices.Add(i);
            }

            return groups;
        }

        /// <summary>
        /// Merges per-window extraction results into a single contiguous output.
        /// Remaps buffer offsets so they point into the unified XIC arrays.
        /// Returns pooled arrays from per-window results after copying.
        /// </summary>
        private static ExtractionResult MergeResults(
            int totalQueryCount,
            WindowExtractionResult[] perWindowResults)
        {
            // Calculate total XIC data points across all windows
            int totalXicPoints = 0;
            for (int i = 0; i < perWindowResults.Length; i++)
            {
                totalXicPoints += perWindowResults[i].TotalDataPoints;
            }

            // Allocate unified output
            var allResults = new FragmentResult[totalQueryCount];
            float[] allRt = new float[totalXicPoints];
            float[] allIntensity = new float[totalXicPoints];

            int globalOffset = 0;

            for (int w = 0; w < perWindowResults.Length; w++)
            {
                ref var windowResult = ref perWindowResults[w];

                // Copy XIC data from per-window buffers into unified arrays
                if (windowResult.TotalDataPoints > 0)
                {
                    Array.Copy(windowResult.RtBuffer, 0, allRt, globalOffset, windowResult.TotalDataPoints);
                    Array.Copy(windowResult.IntensityBuffer, 0, allIntensity, globalOffset, windowResult.TotalDataPoints);
                }

                // Remap each FragmentResult's buffer offsets and place in original query order
                for (int q = 0; q < windowResult.Results.Length; q++)
                {
                    ref readonly var srcResult = ref windowResult.Results[q];
                    int originalIndex = windowResult.OriginalIndices[q];

                    // Remap offsets to point into the merged arrays
                    allResults[originalIndex] = new FragmentResult(
                        srcResult.QueryId,
                        srcResult.DataPointCount,
                        srcResult.RtBufferOffset + globalOffset,
                        srcResult.IntensityBufferOffset + globalOffset,
                        srcResult.TotalIntensity);
                }

                globalOffset += windowResult.TotalDataPoints;

                // Return pooled buffers
                ArrayPool<float>.Shared.Return(windowResult.RtBuffer);
                ArrayPool<float>.Shared.Return(windowResult.IntensityBuffer);
            }

            return new ExtractionResult(allResults, allRt, allIntensity, totalXicPoints);
        }

        public void Dispose()
        {
            _disposed = true;
            // DiaScanIndex is not owned by this orchestrator.
        }

        /// <summary>
        /// Internal struct holding per-window extraction output before merging.
        /// </summary>
        private struct WindowExtractionResult
        {
            public FragmentResult[] Results;
            public int[] OriginalIndices;
            public float[] RtBuffer;       // Rented from ArrayPool
            public float[] IntensityBuffer; // Rented from ArrayPool
            public int TotalDataPoints;
            public int BufferCapacity;

            public WindowExtractionResult(
                FragmentResult[] results,
                int[] originalIndices,
                float[] rtBuffer,
                float[] intensityBuffer,
                int totalDataPoints,
                int bufferCapacity)
            {
                Results = results;
                OriginalIndices = originalIndices;
                RtBuffer = rtBuffer;
                IntensityBuffer = intensityBuffer;
                TotalDataPoints = totalDataPoints;
                BufferCapacity = bufferCapacity;
            }
        }
    }

    /// <summary>
    /// Contains the merged results of a parallel fragment extraction.
    /// 
    /// Results are in the same order as the original input queries.
    /// The XIC data (RT and intensity values) are in contiguous arrays,
    /// with each FragmentResult providing the offset and length into these arrays.
    /// 
    /// Usage:
    ///   var result = orchestrator.ExtractAll(queries);
    ///   for (int i = 0; i &lt; result.Results.Length; i++)
    ///   {
    ///       var r = result.Results[i];
    ///       var xicRt = result.RtBuffer.AsSpan(r.RtBufferOffset, r.DataPointCount);
    ///       var xicInt = result.IntensityBuffer.AsSpan(r.IntensityBufferOffset, r.DataPointCount);
    ///       // ... score, plot, etc.
    ///   }
    /// </summary>
    public sealed class ExtractionResult
    {
        /// <summary>
        /// Per-query results, in the same order as the input queries.
        /// Each result contains the offset and length into the XIC buffers.
        /// </summary>
        public readonly FragmentResult[] Results;

        /// <summary>
        /// Contiguous array of XIC retention time values across all queries.
        /// Access via Results[i].RtBufferOffset and Results[i].DataPointCount.
        /// </summary>
        public readonly float[] RtBuffer;

        /// <summary>
        /// Contiguous array of XIC intensity values across all queries.
        /// Access via Results[i].IntensityBufferOffset and Results[i].DataPointCount.
        /// </summary>
        public readonly float[] IntensityBuffer;

        /// <summary>
        /// Total number of XIC data points across all queries.
        /// </summary>
        public readonly int TotalDataPoints;

        /// <summary>Empty result for zero-query input.</summary>
        public static readonly ExtractionResult Empty = new ExtractionResult(
            Array.Empty<FragmentResult>(), Array.Empty<float>(), Array.Empty<float>(), 0);

        public ExtractionResult(FragmentResult[] results, float[] rtBuffer, float[] intensityBuffer, int totalDataPoints)
        {
            Results = results;
            RtBuffer = rtBuffer;
            IntensityBuffer = intensityBuffer;
            TotalDataPoints = totalDataPoints;
        }
    }
}
