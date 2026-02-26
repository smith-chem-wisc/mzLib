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
        /// Used as fallback when iRT calibration is unavailable or unreliable.
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

        // ── iRT Calibration Parameters ──────────────────────────────────────

        /// <summary>
        /// Whether to use iRT-based calibration for RT windowing and scoring.
        /// When false, falls back to fixed RtToleranceMinutes windows.
        /// Default: true.
        /// </summary>
        public bool UseIrtCalibration { get; set; } = true;

        /// <summary>
        /// Initial broad window half-width in iRT units for Phase 1 (uncalibrated) search.
        /// Must be large enough to capture true matches despite unknown calibration.
        /// Default: 20 iRT units.
        /// </summary>
        public double InitialIrtWindow { get; set; } = 20.0;

        /// <summary>
        /// Number of standard deviations (k) for the calibrated RT window.
        /// Window = k * σ_iRT in iRT space.
        /// Default: 3.0.
        /// </summary>
        public double CalibratedWindowSigmaMultiplier { get; set; } = 3.0;

        /// <summary>
        /// Minimum spectral score for a match to be used as an anchor in calibration fitting.
        /// Higher values produce cleaner calibration at the cost of fewer anchors.
        /// Default: 0.5 (dot product).
        /// </summary>
        public float CalibrationAnchorMinScore { get; set; } = 0.5f;

        /// <summary>
        /// Weight (λ) for the RT score term in the combined score:
        /// finalScore = spectralScore + λ * rtScore
        /// 
        /// λ = 0 disables RT scoring. Typical values: 0.1–1.0.
        /// Default: 0.5.
        /// </summary>
        public double RtScoreLambda { get; set; } = 0.5;

        /// <summary>
        /// Maximum number of calibration refinement iterations.
        /// Each iteration: re-select anchors → refit → re-window → re-score.
        /// Default: 3.
        /// </summary>
        public int MaxCalibrationIterations { get; set; } = 3;

        /// <summary>
        /// Convergence threshold for slope (|a_new - a_old| &lt; epsilon stops iteration).
        /// Default: 0.001.
        /// </summary>
        public double CalibrationConvergenceEpsilon { get; set; } = 0.001;

        /// <summary>
        /// Resolves effective thread count, replacing -1 with processor count.
        /// </summary>
        public int EffectiveMaxThreads =>
            MaxThreads <= 0 ? System.Environment.ProcessorCount : MaxThreads;
    }
}
