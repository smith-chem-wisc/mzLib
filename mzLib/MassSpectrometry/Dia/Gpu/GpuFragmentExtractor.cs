// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// ╔══════════════════════════════════════════════════════════════════════╗
// ║  ILGPU DEPENDENCY REQUIRED                                         ║
// ║                                                                    ║
// ║  This file requires:                                               ║
// ║    dotnet add package ILGPU                                        ║
// ║    dotnet add package ILGPU.Algorithms                             ║
// ║                                                                    ║
// ║  Without these packages, this ONE file will not compile.           ║
// ║  All other DIA engine files compile and run without ILGPU.         ║
// ║  CpuFragmentExtractor is always available as the fallback.         ║
// ║                                                                    ║
// ║  Cross-platform notes:                                             ║
// ║    Windows + NVIDIA GPU   → CUDA backend (best performance)        ║
// ║    Linux   + NVIDIA GPU   → CUDA backend (best performance)        ║
// ║    Windows/Linux + AMD    → OpenCL backend                         ║
// ║    macOS (any)            → Not supported, use CpuFragmentExtractor║
// ╚══════════════════════════════════════════════════════════════════════╝

using System;
using System.Linq;
using ILGPU;
using ILGPU.Runtime;
using ILGPU.Runtime.Cuda;
using ILGPU.Runtime.OpenCL;

namespace MassSpectrometry.Dia
{
    // ══════════════════════════════════════════════════════════════════════
    //  Kernel parameter structs
    //
    //  .NET Action<> supports at most 16 type parameters.
    //  ILGPU's LoadAutoGroupedStreamKernel inherits this limit.
    //  We pack related ArrayViews into structs to keep the parameter
    //  count at 7 (Index1D + 6 structs), well within the limit.
    //  ILGPU handles passing these structs by value to GPU threads.
    // ══════════════════════════════════════════════════════════════════════

    /// <summary>Contiguous peak data arrays (uploaded once at construction).</summary>
    public readonly struct GpuPeakData
    {
        public readonly ArrayView1D<float, Stride1D.Dense> AllMz;
        public readonly ArrayView1D<float, Stride1D.Dense> AllIntensity;

        public GpuPeakData(
            ArrayView1D<float, Stride1D.Dense> allMz,
            ArrayView1D<float, Stride1D.Dense> allIntensity)
        {
            AllMz = allMz;
            AllIntensity = allIntensity;
        }
    }

    /// <summary>Per-scan metadata arrays (uploaded once at construction).</summary>
    public readonly struct GpuScanMeta
    {
        public readonly ArrayView1D<int, Stride1D.Dense> Offsets;
        public readonly ArrayView1D<int, Stride1D.Dense> Lengths;
        public readonly ArrayView1D<float, Stride1D.Dense> Rts;

        public GpuScanMeta(
            ArrayView1D<int, Stride1D.Dense> offsets,
            ArrayView1D<int, Stride1D.Dense> lengths,
            ArrayView1D<float, Stride1D.Dense> rts)
        {
            Offsets = offsets;
            Lengths = lengths;
            Rts = rts;
        }
    }

    /// <summary>Window-to-scan mapping arrays (uploaded once at construction).</summary>
    public readonly struct GpuWindowMap
    {
        public readonly ArrayView1D<int, Stride1D.Dense> ScanStart;
        public readonly ArrayView1D<int, Stride1D.Dense> ScanCount;

        public GpuWindowMap(
            ArrayView1D<int, Stride1D.Dense> scanStart,
            ArrayView1D<int, Stride1D.Dense> scanCount)
        {
            ScanStart = scanStart;
            ScanCount = scanCount;
        }
    }

    /// <summary>Per-query parameters (uploaded per batch).</summary>
    public readonly struct GpuQueryData
    {
        public readonly ArrayView1D<float, Stride1D.Dense> MzLow;
        public readonly ArrayView1D<float, Stride1D.Dense> MzHigh;
        public readonly ArrayView1D<float, Stride1D.Dense> RtMin;
        public readonly ArrayView1D<float, Stride1D.Dense> RtMax;
        public readonly ArrayView1D<int, Stride1D.Dense> WindowId;

        public GpuQueryData(
            ArrayView1D<float, Stride1D.Dense> mzLow,
            ArrayView1D<float, Stride1D.Dense> mzHigh,
            ArrayView1D<float, Stride1D.Dense> rtMin,
            ArrayView1D<float, Stride1D.Dense> rtMax,
            ArrayView1D<int, Stride1D.Dense> windowId)
        {
            MzLow = mzLow;
            MzHigh = mzHigh;
            RtMin = rtMin;
            RtMax = rtMax;
            WindowId = windowId;
        }
    }

    /// <summary>Per-batch output buffers (downloaded after kernel).</summary>
    public readonly struct GpuOutputBuffers
    {
        public readonly ArrayView1D<float, Stride1D.Dense> Intensity;
        public readonly ArrayView1D<float, Stride1D.Dense> Rt;

        public GpuOutputBuffers(
            ArrayView1D<float, Stride1D.Dense> intensity,
            ArrayView1D<float, Stride1D.Dense> rt)
        {
            Intensity = intensity;
            Rt = rt;
        }
    }

    /// <summary>Scalar configuration values for the kernel.</summary>
    public readonly struct GpuKernelConfig
    {
        public readonly int MaxScansPerWindow;
        public readonly int TotalScanCount;

        public GpuKernelConfig(int maxScansPerWindow, int totalScanCount)
        {
            MaxScansPerWindow = maxScansPerWindow;
            TotalScanCount = totalScanCount;
        }
    }

    /// <summary>
    /// GPU-accelerated fragment ion extraction from DIA data via ILGPU.
    /// 
    /// Implements IFragmentExtractor identically to CpuFragmentExtractor, enabling
    /// transparent substitution through FragmentExtractorFactory and the orchestrator.
    /// 
    /// Architecture:
    ///   Construction (one-time):
    ///     - Creates ILGPU Context, selects best Accelerator (CUDA > OpenCL)
    ///     - Uploads entire SoA peak data to GPU memory
    ///     - Compiles the extraction kernel
    ///     - Pre-allocates pooled GPU buffers for query parameters and output
    ///     - PCIe cost: ~50–200 ms for 15M peaks (120 MB)
    /// 
    ///   Per ExtractBatch call:
    ///     - Reuses pre-allocated GPU buffers (grows only if batch exceeds capacity)
    ///     - Uploads query parameters via pre-allocated host arrays (zero GC pressure)
    ///     - Launches kernel: one GPU thread per (query, scan-in-window) pair
    ///     - Downloads into reusable host arrays (zero GC pressure)
    ///     - Compacts sparse output grid into dense XIC buffers on CPU
    /// 
    /// Buffer pooling eliminates the ~1.8s per-batch overhead from GPU memory
    /// allocation/deallocation that dominated the original implementation.
    /// 
    /// Thread safety: NOT thread-safe. One instance per thread.
    /// </summary>
    public sealed class GpuFragmentExtractor : IFragmentExtractor
    {
        private readonly DiaScanIndex _index;
        private readonly Context _context;
        private readonly Accelerator _accelerator;
        private bool _disposed;

        // ── GPU-resident buffers (uploaded once, read-only during extraction) ──
        private readonly MemoryBuffer1D<float, Stride1D.Dense> _gpuAllMz;
        private readonly MemoryBuffer1D<float, Stride1D.Dense> _gpuAllIntensity;
        private readonly MemoryBuffer1D<int, Stride1D.Dense> _gpuScanOffsets;
        private readonly MemoryBuffer1D<int, Stride1D.Dense> _gpuScanLengths;
        private readonly MemoryBuffer1D<float, Stride1D.Dense> _gpuScanRts;
        private readonly MemoryBuffer1D<int, Stride1D.Dense> _gpuWindowScanStart;
        private readonly MemoryBuffer1D<int, Stride1D.Dense> _gpuWindowScanCount;
        private readonly int _maxWindowId;
        private readonly int _maxScansPerWindow;

        // ── Pooled per-batch GPU buffers (allocated once, grown when needed) ──
        // Query parameter buffers (size = current query capacity)
        private MemoryBuffer1D<float, Stride1D.Dense> _poolMzLow;
        private MemoryBuffer1D<float, Stride1D.Dense> _poolMzHigh;
        private MemoryBuffer1D<float, Stride1D.Dense> _poolRtMin;
        private MemoryBuffer1D<float, Stride1D.Dense> _poolRtMax;
        private MemoryBuffer1D<int, Stride1D.Dense> _poolWindowId;
        private int _pooledQueryCapacity;

        // Output buffers (size = current grid capacity = queryCapacity × maxScansPerWindow)
        private MemoryBuffer1D<float, Stride1D.Dense> _poolOutputIntensity;
        private MemoryBuffer1D<float, Stride1D.Dense> _poolOutputRt;
        private int _pooledGridCapacity;

        // ── Reusable host-side arrays (avoid GC allocations per batch) ────────
        private float[] _hostMzLow;
        private float[] _hostMzHigh;
        private float[] _hostRtMin;
        private float[] _hostRtMax;
        private int[] _hostWindowId;
        private float[] _hostOutputIntensity;
        private float[] _hostOutputRt;

        // Compiled kernel (7 params: Index1D + 6 structs — within Action<> limit)
        private readonly Action<Index1D,
            GpuPeakData, GpuScanMeta, GpuWindowMap,
            GpuQueryData, GpuOutputBuffers, GpuKernelConfig
        > _extractKernel;

        /// <summary>Name of the GPU device (for diagnostics/logging).</summary>
        public string DeviceName => _accelerator.Name;

        /// <summary>
        /// Creates a GPU fragment extractor, uploading all peak data to GPU memory.
        /// Pre-allocates pooled buffers for an initial batch size.
        /// Throws InvalidOperationException if no compatible GPU device is found.
        /// </summary>
        /// <param name="index">The DIA scan index to extract from.</param>
        /// <param name="initialQueryCapacity">
        /// Initial capacity for pooled query buffers. Grows automatically if exceeded.
        /// Default 1024 is suitable for most per-window batch sizes.
        /// </param>
        public GpuFragmentExtractor(DiaScanIndex index, int initialQueryCapacity = 1024)
        {
            _index = index ?? throw new ArgumentNullException(nameof(index));

            _context = Context.Create(builder =>
            {
                builder.Cuda().OpenCL().Optimize(OptimizationLevel.O2);
            });

            _accelerator = SelectBestAccelerator();

            // ── Upload SoA peak data to GPU ───────────────────────────────────
            float[] allMzArr = index.AllMz.ToArray();
            float[] allIntArr = index.AllIntensity.ToArray();
            _gpuAllMz = _accelerator.Allocate1D<float>(allMzArr.Length);
            _gpuAllMz.CopyFromCPU(allMzArr);
            _gpuAllIntensity = _accelerator.Allocate1D<float>(allIntArr.Length);
            _gpuAllIntensity.CopyFromCPU(allIntArr);

            int[] scanOffsets = index.AllScanOffsets.ToArray();
            int[] scanLengths = index.AllScanLengths.ToArray();
            float[] scanRts = index.AllScanRts.ToArray();

            _gpuScanOffsets = _accelerator.Allocate1D<int>(scanOffsets.Length);
            _gpuScanOffsets.CopyFromCPU(scanOffsets);
            _gpuScanLengths = _accelerator.Allocate1D<int>(scanLengths.Length);
            _gpuScanLengths.CopyFromCPU(scanLengths);
            _gpuScanRts = _accelerator.Allocate1D<float>(scanRts.Length);
            _gpuScanRts.CopyFromCPU(scanRts);

            // ── Build flat window→scan mapping ────────────────────────────────
            var windowIds = index.GetWindowIds();
            _maxWindowId = -1;
            _maxScansPerWindow = 0;

            foreach (int wid in windowIds)
            {
                if (wid > _maxWindowId) _maxWindowId = wid;
                if (index.TryGetScanRangeForWindow(wid, out _, out int count))
                    if (count > _maxScansPerWindow) _maxScansPerWindow = count;
            }

            int windowArraySize = _maxWindowId + 1;
            int[] windowScanStart = new int[windowArraySize];
            int[] windowScanCount = new int[windowArraySize];

            foreach (int wid in windowIds)
            {
                if (index.TryGetScanRangeForWindow(wid, out int start, out int count))
                {
                    windowScanStart[wid] = start;
                    windowScanCount[wid] = count;
                }
            }

            _gpuWindowScanStart = _accelerator.Allocate1D<int>(windowArraySize);
            _gpuWindowScanStart.CopyFromCPU(windowScanStart);
            _gpuWindowScanCount = _accelerator.Allocate1D<int>(windowArraySize);
            _gpuWindowScanCount.CopyFromCPU(windowScanCount);

            // ── Pre-allocate pooled buffers ───────────────────────────────────
            EnsureQueryCapacity(initialQueryCapacity);

            // ── Compile kernel ────────────────────────────────────────────────
            _extractKernel = _accelerator.LoadAutoGroupedStreamKernel<
                Index1D, GpuPeakData, GpuScanMeta, GpuWindowMap,
                GpuQueryData, GpuOutputBuffers, GpuKernelConfig>(ExtractKernel);

            _accelerator.Synchronize();
        }

        /// <summary>
        /// Ensures pooled GPU buffers and host arrays are large enough for the given
        /// query count. Only reallocates when capacity is exceeded — otherwise a no-op.
        /// Uses 1.5× growth factor to amortize reallocation cost.
        /// </summary>
        private void EnsureQueryCapacity(int queryCount)
        {
            if (queryCount <= _pooledQueryCapacity) return;

            // Grow with 1.5× factor to avoid frequent reallocations
            int newCapacity = Math.Max(queryCount, (int)(_pooledQueryCapacity * 1.5));
            newCapacity = Math.Max(newCapacity, 128); // minimum useful size

            // Dispose old GPU query buffers (if any)
            _poolMzLow?.Dispose();
            _poolMzHigh?.Dispose();
            _poolRtMin?.Dispose();
            _poolRtMax?.Dispose();
            _poolWindowId?.Dispose();

            // Allocate new GPU query buffers
            _poolMzLow = _accelerator.Allocate1D<float>(newCapacity);
            _poolMzHigh = _accelerator.Allocate1D<float>(newCapacity);
            _poolRtMin = _accelerator.Allocate1D<float>(newCapacity);
            _poolRtMax = _accelerator.Allocate1D<float>(newCapacity);
            _poolWindowId = _accelerator.Allocate1D<int>(newCapacity);
            _pooledQueryCapacity = newCapacity;

            // Allocate host-side upload arrays
            _hostMzLow = new float[newCapacity];
            _hostMzHigh = new float[newCapacity];
            _hostRtMin = new float[newCapacity];
            _hostRtMax = new float[newCapacity];
            _hostWindowId = new int[newCapacity];

            // Reallocate output grid buffers to match
            int newGridCapacity = newCapacity * _maxScansPerWindow;
            EnsureGridCapacity(newGridCapacity);
        }

        /// <summary>
        /// Ensures output grid GPU buffers and host download arrays are large enough.
        /// Only reallocates when capacity is exceeded.
        /// </summary>
        private void EnsureGridCapacity(int gridSize)
        {
            if (gridSize <= _pooledGridCapacity) return;

            int newCapacity = Math.Max(gridSize, (int)(_pooledGridCapacity * 1.5));
            newCapacity = Math.Max(newCapacity, 128);

            _poolOutputIntensity?.Dispose();
            _poolOutputRt?.Dispose();

            _poolOutputIntensity = _accelerator.Allocate1D<float>(newCapacity);
            _poolOutputRt = _accelerator.Allocate1D<float>(newCapacity);
            _pooledGridCapacity = newCapacity;

            // Host-side download arrays
            _hostOutputIntensity = new float[newCapacity];
            _hostOutputRt = new float[newCapacity];
        }

        // ══════════════════════════════════════════════════════════════════════
        //  GPU KERNEL
        //
        //  Thread layout: flat 1D, one thread per (query, scan-in-window) pair.
        //  threadIndex = queryIdx * maxScansPerWindow + scanOffsetInWindow
        //
        //  Binary search + intensity accumulation — identical logic to CPU path.
        // ══════════════════════════════════════════════════════════════════════

        private static void ExtractKernel(
            Index1D threadIndex,
            GpuPeakData peaks,
            GpuScanMeta scans,
            GpuWindowMap windows,
            GpuQueryData queries,
            GpuOutputBuffers output,
            GpuKernelConfig config)
        {
            int queryIdx = threadIndex / config.MaxScansPerWindow;
            int scanOffsetInWindow = threadIndex % config.MaxScansPerWindow;

            if (queryIdx >= queries.MzLow.Length) return;

            int wid = queries.WindowId[queryIdx];
            if (wid < 0 || wid >= windows.ScanStart.Length) return;

            int wStart = windows.ScanStart[wid];
            int wCount = windows.ScanCount[wid];
            if (scanOffsetInWindow >= wCount) return;

            int scanIdx = wStart + scanOffsetInWindow;
            if (scanIdx >= config.TotalScanCount) return;

            // RT filtering
            float rt = scans.Rts[scanIdx];
            if (rt < queries.RtMin[queryIdx] || rt > queries.RtMax[queryIdx]) return;

            float mzLow = queries.MzLow[queryIdx];
            float mzHigh = queries.MzHigh[queryIdx];

            int peakStart = scans.Offsets[scanIdx];
            int peakCount = scans.Lengths[scanIdx];
            if (peakCount == 0) return;

            // Binary search: first index where mz >= mzLow
            int lo = 0;
            int hi = peakCount;
            while (lo < hi)
            {
                int mid = lo + ((hi - lo) >> 1);
                if (peaks.AllMz[peakStart + mid] < mzLow)
                    lo = mid + 1;
                else
                    hi = mid;
            }

            // Sum intensities of peaks in [mzLow, mzHigh]
            float intensity = 0f;
            for (int p = lo; p < peakCount; p++)
            {
                float mz = peaks.AllMz[peakStart + p];
                if (mz > mzHigh) break;
                intensity += peaks.AllIntensity[peakStart + p];
            }

            if (intensity > 0f)
            {
                int outputIdx = queryIdx * config.MaxScansPerWindow + scanOffsetInWindow;
                output.Intensity[outputIdx] = intensity;
                output.Rt[outputIdx] = rt;
            }
        }

        // ══════════════════════════════════════════════════════════════════════
        //  IFragmentExtractor implementation
        // ══════════════════════════════════════════════════════════════════════

        /// <inheritdoc/>
        public int ExtractBatch(
            ReadOnlySpan<FragmentQuery> queries,
            Span<FragmentResult> results,
            Span<float> rtBuffer,
            Span<float> intensityBuffer)
        {
            if (queries.Length != results.Length)
                throw new ArgumentException("Queries and results spans must have the same length.");

            int queryCount = queries.Length;
            if (queryCount == 0) return 0;

            // ── Ensure pooled buffers are large enough ────────────────────────
            EnsureQueryCapacity(queryCount);
            int gridSize = queryCount * _maxScansPerWindow;
            EnsureGridCapacity(gridSize);

            // ── Fill reusable host arrays ─────────────────────────────────────
            for (int i = 0; i < queryCount; i++)
            {
                float tol = queries[i].TolerancePpm / 1_000_000f;
                _hostMzLow[i] = queries[i].TargetMz * (1f - tol);
                _hostMzHigh[i] = queries[i].TargetMz * (1f + tol);
                _hostRtMin[i] = queries[i].RtMin;
                _hostRtMax[i] = queries[i].RtMax;
                _hostWindowId[i] = queries[i].WindowId;
            }

            // ── Upload query data into pooled GPU buffers ─────────────────────
            // Host arrays and GPU buffers are both sized to _pooledQueryCapacity.
            // Upload full arrays, then pass SubViews to the kernel for bounds checking.
            _poolMzLow.CopyFromCPU(_hostMzLow);
            _poolMzHigh.CopyFromCPU(_hostMzHigh);
            _poolRtMin.CopyFromCPU(_hostRtMin);
            _poolRtMax.CopyFromCPU(_hostRtMax);
            _poolWindowId.CopyFromCPU(_hostWindowId);

            // Zero only the portion of the output grid we'll use
            _poolOutputIntensity.View.SubView(0, gridSize).MemSetToZero();
            _poolOutputRt.View.SubView(0, gridSize).MemSetToZero();

            // ── Launch kernel using pooled buffer views ────────────────────────
            // Pass SubViews sized to actual queryCount so the kernel bounds-checks work
            var peakViews = new GpuPeakData(_gpuAllMz.View, _gpuAllIntensity.View);
            var scanViews = new GpuScanMeta(_gpuScanOffsets.View, _gpuScanLengths.View, _gpuScanRts.View);
            var windowViews = new GpuWindowMap(_gpuWindowScanStart.View, _gpuWindowScanCount.View);
            var queryViews = new GpuQueryData(
                _poolMzLow.View.SubView(0, queryCount),
                _poolMzHigh.View.SubView(0, queryCount),
                _poolRtMin.View.SubView(0, queryCount),
                _poolRtMax.View.SubView(0, queryCount),
                _poolWindowId.View.SubView(0, queryCount));
            var outputViews = new GpuOutputBuffers(
                _poolOutputIntensity.View.SubView(0, gridSize),
                _poolOutputRt.View.SubView(0, gridSize));
            var config = new GpuKernelConfig(_maxScansPerWindow, _index.ScanCount);

            _extractKernel(gridSize, peakViews, scanViews, windowViews, queryViews, outputViews, config);
            _accelerator.Synchronize();

            // ── Download into reusable host arrays ────────────────────────────
            _poolOutputIntensity.CopyToCPU(_hostOutputIntensity);
            _poolOutputRt.CopyToCPU(_hostOutputRt);

            // ── Compact sparse grid into dense output buffers ─────────────────
            int totalDataPoints = 0;
            for (int q = 0; q < queryCount; q++)
            {
                int dataPointStart = totalDataPoints;
                float totalIntensity = 0f;
                int gridBase = q * _maxScansPerWindow;

                for (int s = 0; s < _maxScansPerWindow; s++)
                {
                    float inten = _hostOutputIntensity[gridBase + s];
                    if (inten > 0f)
                    {
                        rtBuffer[totalDataPoints] = _hostOutputRt[gridBase + s];
                        intensityBuffer[totalDataPoints] = inten;
                        totalDataPoints++;
                        totalIntensity += inten;
                    }
                }

                int dataPointCount = totalDataPoints - dataPointStart;
                results[q] = new FragmentResult(
                    queries[q].QueryId, dataPointCount,
                    dataPointStart, dataPointStart, totalIntensity);
            }

            return totalDataPoints;
        }

        private Accelerator SelectBestAccelerator()
        {
            foreach (var device in _context.Devices)
                if (device.AcceleratorType == AcceleratorType.Cuda)
                    return device.CreateAccelerator(_context);

            foreach (var device in _context.Devices)
                if (device.AcceleratorType == AcceleratorType.OpenCL)
                    return device.CreateAccelerator(_context);

            throw new InvalidOperationException(
                "No compatible GPU accelerator found. Ensure NVIDIA drivers (CUDA) or " +
                "OpenCL runtime is installed. Use CpuFragmentExtractor as fallback.");
        }

        public void Dispose()
        {
            if (!_disposed)
            {
                _disposed = true;

                // Pooled per-batch buffers
                _poolMzLow?.Dispose();
                _poolMzHigh?.Dispose();
                _poolRtMin?.Dispose();
                _poolRtMax?.Dispose();
                _poolWindowId?.Dispose();
                _poolOutputIntensity?.Dispose();
                _poolOutputRt?.Dispose();

                // Static data buffers
                _gpuAllMz?.Dispose();
                _gpuAllIntensity?.Dispose();
                _gpuScanOffsets?.Dispose();
                _gpuScanLengths?.Dispose();
                _gpuScanRts?.Dispose();
                _gpuWindowScanStart?.Dispose();
                _gpuWindowScanCount?.Dispose();

                _accelerator?.Dispose();
                _context?.Dispose();
            }
        }
    }
}