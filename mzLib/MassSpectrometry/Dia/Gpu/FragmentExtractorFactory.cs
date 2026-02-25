// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Creates the best available IFragmentExtractor for a given DiaScanIndex.
    /// 
    /// Automatically selects GPU when available, falls back to CPU otherwise.
    /// Uses reflection to instantiate GpuFragmentExtractor so this class (and the
    /// rest of the core mzLib assembly) has no compile-time dependency on ILGPU.
    /// 
    /// Usage (single extractor):
    ///   using var extractor = FragmentExtractorFactory.Create(index);
    /// 
    /// Usage (with orchestrator — one extractor per thread):
    ///   var orch = new DiaExtractionOrchestrator(index,
    ///       FragmentExtractorFactory.CreateFactory());
    /// </summary>
    public static class FragmentExtractorFactory
    {
        /// <summary>
        /// Creates the best available fragment extractor.
        /// 
        /// Selection order:
        ///   1. GPU (CUDA or OpenCL) if ILGPU installed + compatible device found
        ///   2. CPU (always available)
        /// 
        /// If GPU creation fails at runtime, falls back to CPU with a diagnostic
        /// message on Console.Error.
        /// </summary>
        /// <param name="index">The DIA scan index to extract from.</param>
        /// <param name="preferCpu">Force CPU even when GPU is available (for benchmarking).</param>
        public static IFragmentExtractor Create(DiaScanIndex index, bool preferCpu = false)
        {
            if (index == null) throw new ArgumentNullException(nameof(index));

            if (!preferCpu && GpuDeviceDetector.IsGpuAvailable)
            {
                try
                {
                    return InstantiateGpuExtractor(index);
                }
                catch (Exception ex)
                {
                    Console.Error.WriteLine(
                        $"[DIA] GPU extractor creation failed, falling back to CPU: {ex.Message}");
                }
            }

            return new CpuFragmentExtractor(index);
        }

        /// <summary>
        /// Returns a factory function for the orchestrator's per-thread extractor creation.
        /// The backend decision is made once; the returned function always creates that type.
        /// </summary>
        /// <param name="preferCpu">Force CPU even when GPU is available.</param>
        public static Func<DiaScanIndex, IFragmentExtractor> CreateFactory(bool preferCpu = false)
        {
            bool useGpu = !preferCpu && GpuDeviceDetector.IsGpuAvailable;

            if (useGpu)
            {
                return idx =>
                {
                    try { return InstantiateGpuExtractor(idx); }
                    catch { return new CpuFragmentExtractor(idx); }
                };
            }

            return idx => new CpuFragmentExtractor(idx);
        }

        /// <summary>
        /// Describes which backend would be selected, without actually creating an extractor.
        /// Useful for logging at startup.
        /// </summary>
        public static string DescribeBackend(bool preferCpu = false)
        {
            if (preferCpu) return "CPU (forced by preference)";
            if (!GpuDeviceDetector.IsGpuAvailable)
                return $"CPU ({GpuDeviceDetector.Description})";
            return $"GPU ({GpuDeviceDetector.DetectedBackend}: {GpuDeviceDetector.Description})";
        }

        /// <summary>
        /// Instantiates GpuFragmentExtractor via reflection.
        /// This keeps the core assembly free of any compile-time reference to ILGPU.
        /// </summary>
        private static IFragmentExtractor InstantiateGpuExtractor(DiaScanIndex index)
        {
            // Look in all loaded assemblies for the type
            var gpuType = Type.GetType(
                "MassSpectrometry.Dia.GpuFragmentExtractor, MassSpectrometry");

            // If not found by qualified name, try scanning loaded assemblies
            if (gpuType == null)
            {
                foreach (var asm in AppDomain.CurrentDomain.GetAssemblies())
                {
                    gpuType = asm.GetType("MassSpectrometry.Dia.GpuFragmentExtractor");
                    if (gpuType != null) break;
                }
            }

            if (gpuType == null)
            {
                throw new InvalidOperationException(
                    "GpuFragmentExtractor type not found. Ensure ILGPU NuGet packages are " +
                    "installed and GpuFragmentExtractor.cs is compiled into the project.");
            }

            var instance = Activator.CreateInstance(gpuType, index);
            if (instance is IFragmentExtractor extractor)
                return extractor;

            throw new InvalidOperationException(
                "GpuFragmentExtractor was found but does not implement IFragmentExtractor.");
        }
    }
}
