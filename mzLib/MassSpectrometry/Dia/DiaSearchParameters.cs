// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Configuration parameters for DIA fragment extraction and scoring.
    /// Passed through the entire pipeline: query generation → extraction → scoring → result assembly.
    /// </summary>
    public class DiaSearchParameters
    {
        /// <summary>
        /// Fragment m/z tolerance in parts-per-million for XIC extraction.
        /// Default: 20 ppm (typical for Orbitrap MS2).
        /// </summary>
        public float PpmTolerance { get; set; } = 20f;

        /// <summary>
        /// Half-width of the RT window around the predicted/library retention time (in minutes).
        /// Queries are restricted to [RT - tolerance, RT + tolerance].
        /// Used by the uncalibrated Generate() path; ignored when GenerateCalibrated() is used.
        /// Default: 5.0 minutes.
        /// </summary>
        public float RtToleranceMinutes { get; set; } = 5.0f;

        /// <summary>
        /// Minimum number of fragment ions that must yield XIC data points
        /// for a precursor to be reported as a match.
        /// Default: 3.
        /// </summary>
        public int MinFragmentsRequired { get; set; } = 3;

        /// <summary>
        /// Minimum dot product or spectral angle score threshold.
        /// Precursors scoring below this are excluded from results.
        /// Default: 0.0 (no filtering; leave to downstream FDR).
        /// </summary>
        public float MinScoreThreshold { get; set; } = 0.0f;

        /// <summary>
        /// Maximum threads for window-level parallelism.
        /// -1 = use Environment.ProcessorCount.
        /// Default: -1.
        /// </summary>
        public int MaxThreads { get; set; } = -1;

        /// <summary>
        /// Whether to prefer GPU extraction when a compatible device is detected.
        /// Falls back to CPU if no GPU is available regardless of this setting.
        /// Default: false.
        /// </summary>
        public bool PreferGpu { get; set; } = false;

        /// <summary>
        /// Multiplier for the calibration model's sigma (residual standard deviation)
        /// to determine the RT extraction window half-width when using GenerateCalibrated().
        /// 
        /// Window half-width = σ × CalibratedWindowSigmaMultiplier
        /// 
        /// For example, with σ = 0.3 min and multiplier = 3.0:
        ///   window = ±0.9 min (covers 99.7% of a Gaussian distribution)
        /// 
        /// Default: 3.0 (3σ).
        /// </summary>
        public double CalibratedWindowSigmaMultiplier { get; set; } = 3.0;

        /// <summary>
        /// Resolves effective thread count, replacing -1 with processor count.
        /// </summary>
        public int EffectiveMaxThreads =>
            MaxThreads <= 0 ? System.Environment.ProcessorCount : MaxThreads;
    }
}