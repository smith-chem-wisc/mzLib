// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Detects GPU compute availability at runtime without requiring a compile-time
    /// dependency on any GPU library.
    /// 
    /// Cross-platform behavior:
    ///   Windows/Linux + NVIDIA GPU + ILGPU NuGet installed → IsGpuAvailable = true (CUDA)
    ///   Windows/Linux + AMD/Intel GPU + ILGPU NuGet installed → IsGpuAvailable = true (OpenCL)
    ///   macOS (any hardware) → IsGpuAvailable = false (CUDA unsupported, OpenCL deprecated)
    ///   Any platform without ILGPU NuGet → IsGpuAvailable = false
    ///   Any platform without a compatible GPU → IsGpuAvailable = false
    /// 
    /// All probing is via Assembly.Load + reflection. This class compiles and runs
    /// on every .NET platform regardless of whether ILGPU is referenced.
    /// 
    /// Detection runs once (lazily) and the result is cached for the process lifetime.
    /// </summary>
    public static class GpuDeviceDetector
    {
        private static readonly Lazy<GpuDetectionResult> _cached =
            new Lazy<GpuDetectionResult>(Detect);

        /// <summary>True if a usable GPU accelerator was found. Never throws.</summary>
        public static bool IsGpuAvailable => _cached.Value.Available;

        /// <summary>Human-readable description of the GPU or reason for unavailability.</summary>
        public static string Description => _cached.Value.Description;

        /// <summary>Which GPU backend was detected (None, Cuda, or OpenCL).</summary>
        public static GpuBackend DetectedBackend => _cached.Value.Backend;

        private static GpuDetectionResult Detect()
        {
            try
            {
                // On macOS, bail early — no viable GPU compute path in .NET.
                if (System.Runtime.InteropServices.RuntimeInformation.IsOSPlatform(
                    System.Runtime.InteropServices.OSPlatform.OSX))
                {
                    return new GpuDetectionResult(false, GpuBackend.None,
                        "macOS detected. GPU compute not supported (CUDA unavailable, OpenCL deprecated).");
                }

                // Try to load the ILGPU assembly
                var ilgpuAssembly = TryLoadAssembly("ILGPU");
                if (ilgpuAssembly == null)
                {
                    return new GpuDetectionResult(false, GpuBackend.None,
                        "ILGPU not installed. Add the ILGPU NuGet package for optional GPU acceleration.");
                }

                // ILGPU is present. In ILGPU 1.x, CUDA and OpenCL runtimes are
                // bundled inside the main ILGPU assembly (not separate assemblies).
                // Check for the runtime types directly.
                bool hasCuda = ilgpuAssembly.GetType("ILGPU.Runtime.Cuda.CudaAccelerator") != null;
                bool hasOpenCl = ilgpuAssembly.GetType("ILGPU.Runtime.OpenCL.CLAccelerator") != null;

                if (!hasCuda && !hasOpenCl)
                {
                    // Fallback: try separate assemblies (older ILGPU versions)
                    hasCuda = TryLoadAssembly("ILGPU.Runtime.Cuda") != null;
                    hasOpenCl = TryLoadAssembly("ILGPU.Runtime.OpenCL") != null;
                }

                if (!hasCuda && !hasOpenCl)
                {
                    return new GpuDetectionResult(false, GpuBackend.None,
                        "ILGPU found but no CUDA or OpenCL runtime types detected.");
                }

                var backend = hasCuda ? GpuBackend.Cuda : GpuBackend.OpenCL;
                string desc = hasCuda
                    ? "ILGPU with CUDA runtime detected. GPU device selected at extractor creation."
                    : "ILGPU with OpenCL runtime detected. GPU device selected at extractor creation.";

                return new GpuDetectionResult(true, backend, desc);
            }
            catch (Exception ex)
            {
                return new GpuDetectionResult(false, GpuBackend.None,
                    $"GPU detection failed: {ex.GetType().Name}: {ex.Message}");
            }
        }

        private static System.Reflection.Assembly TryLoadAssembly(string name)
        {
            try { return System.Reflection.Assembly.Load(name); }
            catch { return null; }
        }

        internal readonly struct GpuDetectionResult
        {
            public readonly bool Available;
            public readonly GpuBackend Backend;
            public readonly string Description;

            public GpuDetectionResult(bool available, GpuBackend backend, string description)
            {
                Available = available;
                Backend = backend;
                Description = description;
            }
        }
    }

    /// <summary>Which GPU compute backend was detected at runtime.</summary>
    public enum GpuBackend
    {
        /// <summary>No GPU. CPU fallback will be used.</summary>
        None = 0,
        /// <summary>NVIDIA CUDA via ILGPU (best performance, requires NVIDIA GPU + driver).</summary>
        Cuda = 1,
        /// <summary>OpenCL via ILGPU (broader hardware support, slightly lower performance).</summary>
        OpenCL = 2,
    }
}