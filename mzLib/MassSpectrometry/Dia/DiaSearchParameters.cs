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
        /// When true, scoring uses the maximum single-scan intensity per fragment (apex)
        /// instead of the sum of all intensities across scans (TotalIntensity).
        /// 
        /// Apex intensity better represents the fragment's relative abundance at the 
        /// chromatographic peak and correlates more strongly with library relative intensities.
        /// TotalIntensity sums signal over the entire RT window, accumulating noise and
        /// interference from co-eluting peptides.
        /// 
        /// Requires XIC buffers to be passed to AssembleResults. If buffers are not available,
        /// falls back to TotalIntensity regardless of this setting.
        /// 
        /// Default: true.
        /// </summary>
        public bool UseApexIntensityForScoring { get; set; } = true;

        /// <summary>
        /// Resolves effective thread count, replacing -1 with processor count.
        /// </summary>
        public int EffectiveMaxThreads =>
            MaxThreads <= 0 ? System.Environment.ProcessorCount : MaxThreads;
    }
}