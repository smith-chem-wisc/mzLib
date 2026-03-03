// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Abstraction layer for fragment ion extraction from DIA data.
    /// 
    /// This interface decouples the DIA processing orchestration from the compute
    /// implementation. Two implementations are planned:
    ///   - CpuFragmentExtractor: uses binary search + SIMD on CPU
    ///   - GpuFragmentExtractor: batches work onto GPU for parallel extraction
    /// 
    /// Both implementations operate on a DiaScanIndex (SoA layout) and process
    /// batches of FragmentQuery structs, writing results into shared buffers.
    /// 
    /// The caller is responsible for:
    ///   1. Providing the DiaScanIndex to the implementation at construction
    ///   2. Allocating result buffers of sufficient size
    ///   3. Not calling ExtractBatch concurrently on the same instance
    ///      (use one instance per thread, or use window-level parallelism externally)
    /// </summary>
    public interface IFragmentExtractor : IDisposable
    {
        /// <summary>
        /// Extracts fragment ion chromatograms for a batch of queries.
        /// 
        /// For each query, the implementation:
        ///   1. Looks up the scan range for the query's WindowId
        ///   2. Restricts to scans within [RtMin, RtMax]
        ///   3. For each qualifying scan, performs m/z lookup within TargetMz ± TolerancePpm
        ///   4. Writes XIC data points into the result buffers
        ///   5. Populates the corresponding FragmentResult with offsets and summary stats
        /// 
        /// All XIC data points (RT, intensity pairs) are written into rtBuffer and 
        /// intensityBuffer starting at the current write position. The caller must ensure
        /// these buffers are large enough. A conservative estimate is:
        ///   queries.Length * maxScansPerWindow
        /// </summary>
        /// <param name="queries">Batch of fragment queries to extract.</param>
        /// <param name="results">
        /// Output span, same length as queries. Each element describes where the 
        /// extracted XIC data lives in the output buffers.
        /// </param>
        /// <param name="rtBuffer">
        /// Shared buffer for XIC retention time values. Written sequentially.
        /// </param>
        /// <param name="intensityBuffer">
        /// Shared buffer for XIC intensity values. Written sequentially, aligned with rtBuffer.
        /// </param>
        /// <returns>
        /// Total number of XIC data points written across all queries.
        /// This is the high-water mark in rtBuffer/intensityBuffer.
        /// </returns>
        int ExtractBatch(
            ReadOnlySpan<FragmentQuery> queries,
            Span<FragmentResult> results,
            Span<float> rtBuffer,
            Span<float> intensityBuffer);
    }
}
