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
        /// Multiplier for the calibrated RT window width, expressed in units of
        /// the calibration model's residual standard deviation (sigma).
        /// 
        /// After iRT calibration fits a linear model (iRT → experimental RT),
        /// the residual sigma describes the typical RT prediction error.
        /// The calibrated extraction window is: predicted_RT ± (sigma × this multiplier).
        /// 
        /// Default: 3.0 (covers ~99.7% of well-calibrated peptides).
        /// Lower values (e.g., 2.0) give narrower windows and faster extraction
        /// but may miss peptides with unusual RT behavior.
        /// </summary>
        public double CalibratedWindowSigmaMultiplier { get; set; } = 3.0;

        /// <summary>
        /// Controls how extracted fragment intensities are aggregated and scored.
        /// 
        /// Summed: original behavior — sum all XIC intensity per fragment, compare to library.
        /// ConsensusApex: score at the chromatographic apex (best single time point).
        /// TemporalCosine: intensity-weighted average cosine across all RT points.
        /// WeightedTemporalCosineWithTransform: DIA-NN-style with sqrt weighting + cos^N transform.
        /// 
        /// Default: TemporalCosine (good balance of quality and interpretability).
        /// </summary>
        public ScoringStrategy ScoringStrategy { get; set; } = ScoringStrategy.TemporalCosine;

        /// <summary>
        /// Exponent for the nonlinear transform in WeightedTemporalCosineWithTransform mode.
        /// Higher values accentuate good matches and suppress mediocre ones.
        /// 
        /// 1.0 = no transform (identical to TemporalCosine).
        /// 3.0 = DIA-NN default. Improves target-decoy separation.
        /// 
        /// Only used when ScoringStrategy is WeightedTemporalCosineWithTransform.
        /// Default: 3.0.
        /// </summary>
        public float NonlinearPower { get; set; } = 3.0f;

        /// <summary>
        /// Resolves effective thread count, replacing -1 with processor count.
        /// </summary>
        public int EffectiveMaxThreads =>
            MaxThreads <= 0 ? System.Environment.ProcessorCount : MaxThreads;
    }
}