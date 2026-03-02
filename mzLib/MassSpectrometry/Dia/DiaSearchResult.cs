// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;

namespace MassSpectrometry.Dia
{

    /// <summary>
    /// Represents a single DIA precursor identification result with multi-feature scoring.
    /// 
    /// Phase 13 expansion: from 17 features to ~30 features for improved target/decoy
    /// discrimination. Features are organized into categories matching DIA-NN's architecture:
    /// 
    ///   Category 1: Spectral Similarity (cosine, spectral angle, temporal cosine)
    ///   Category 2: Fragment Co-elution (correlations, best-fragment reference curve)
    ///   Category 3: Signal/Noise and Intensity (log intensity, S/N, signal ratio deviation)
    ///   Category 4: Peak Shape (boundary ratio, apex-to-mean, candidate count)
    ///   Category 5: Retention Time Quality (RT deviation, normalized)
    /// 
    /// Lives in MassSpectrometry.Dia (not Omics) because it is tightly coupled to the
    /// DIA engine's SoA extraction output format.
    /// </summary>
    public class DiaSearchResult
    {
        #region Identification

        /// <summary>Peptide sequence with modifications, as stored in LibrarySpectrum</summary>
        public string Sequence { get; }

        /// <summary>Precursor charge state</summary>
        public int ChargeState { get; }

        /// <summary>Precursor m/z (from library)</summary>
        public double PrecursorMz { get; }

        /// <summary>DIA isolation window ID this precursor was extracted from</summary>
        public int WindowId { get; }

        /// <summary>Whether this is a decoy identification (from decoy library)</summary>
        public bool IsDecoy { get; }

        #endregion

        #region Category 1: Spectral Similarity Scores

        /// <summary>
        /// Normalized dot product between library and extracted fragment intensities.
        /// Range [0, 1], higher is better. NaN if insufficient fragments.
        /// </summary>
        public float DotProductScore { get; set; }

        /// <summary>
        /// Spectral angle score: 1 - (2/π) × arccos(normalized dot product).
        /// Range [0, 1], higher is better. More sensitive than DP in high-similarity range.
        /// </summary>
        public float SpectralAngleScore { get; set; }

        /// <summary>
        /// Raw cosine value before clamping/transform. Used for downstream spectral angle.
        /// </summary>
        public float RawCosine { get; set; }

        /// <summary>
        /// Cosine similarity at the apex time point (full-window apex).
        /// </summary>
        public float ApexDotProductScore { get; set; }

        /// <summary>
        /// Intensity-weighted average cosine across all time points in the window.
        /// More robust than single-scan apex scoring.
        /// </summary>
        public float TemporalCosineScore { get; set; }

        #endregion

        #region Category 2: Fragment Co-elution / Correlation

        /// <summary>Mean pairwise Pearson correlation between fragment XICs (full window).</summary>
        public float MeanFragmentCorrelation { get; set; }

        /// <summary>Minimum pairwise Pearson correlation between fragment XICs (full window).</summary>
        public float MinFragmentCorrelation { get; set; }

        /// <summary>Mean pairwise Pearson correlation within detected peak boundaries.</summary>
        public float PeakMeanFragCorrelation { get; set; }

        /// <summary>Minimum pairwise Pearson correlation within detected peak boundaries.</summary>
        public float PeakMinFragCorrelation { get; set; }

        // ── Phase 13: Best-Fragment Reference Curve (DIA-NN core innovation) ──

        /// <summary>
        /// Sum of Pearson correlations between the best fragment and all other fragments.
        /// The "best fragment" is the one least affected by interference (max correlation sum).
        /// Higher values indicate cleaner elution with less interference.
        /// </summary>
        public float BestFragCorrelationSum { get; set; }

        /// <summary>
        /// Median correlation of individual fragments with the best reference fragment.
        /// True peptides: high (all fragments track reference). Decoys: low (random).
        /// </summary>
        public float MedianFragRefCorr { get; set; }

        /// <summary>
        /// Minimum correlation of any fragment with the best reference fragment.
        /// Catches the single worst-interfered fragment.
        /// </summary>
        public float MinFragRefCorr { get; set; }

        /// <summary>
        /// Standard deviation of fragment-vs-reference correlations.
        /// True peptides: low (all fragments correlate similarly). Interference: high (some fragments diverge).
        /// </summary>
        public float StdFragRefCorr { get; set; }

        /// <summary>
        /// Cosine similarity at apex, weighted by each fragment's correlation with the best reference.
        /// Down-weights interfered fragments, boosting score for true peptides.
        /// </summary>
        public float BestFragWeightedCosine { get; set; }

        /// <summary>
        /// Index of the best (least-interfered) fragment. For diagnostic purposes.
        /// </summary>
        public int BestFragIndex { get; set; }

        // ── Phase 13: Smoothed Correlations ──

        /// <summary>
        /// Mean pairwise Pearson correlation on min-of-3-consecutive smoothed XICs.
        /// Reduces noise spike artifacts that can inflate/deflate raw correlations.
        /// </summary>
        public float SmoothedMeanFragCorr { get; set; }

        /// <summary>
        /// Minimum pairwise correlation on smoothed XICs.
        /// </summary>
        public float SmoothedMinFragCorr { get; set; }

        #endregion

        #region Category 3: Signal/Noise and Intensity

        /// <summary>
        /// Log2 of total extracted intensity within the peak/window.
        /// Raw signal strength is independently discriminative.
        /// </summary>
        public float LogTotalIntensity { get; set; }

        /// <summary>
        /// Log2(signal / noise) where signal = apex total intensity,
        /// noise = median non-peak total intensity.
        /// </summary>
        public float Log2SignalToNoise { get; set; }

        // ── Phase 13: Per-Fragment Signal Ratio Deviation ──

        /// <summary>
        /// Mean absolute log2-ratio deviation: |log2(obs_fraction / lib_fraction)| per fragment.
        /// True peptides: ~0 (observed ratios match library). Interference: large (specific fragments inflated).
        /// LOWER is better for targets.
        /// </summary>
        public float MeanSignalRatioDeviation { get; set; }

        /// <summary>
        /// Maximum absolute log2-ratio deviation across fragments.
        /// Identifies the single most interfered fragment.
        /// </summary>
        public float MaxSignalRatioDeviation { get; set; }

        /// <summary>
        /// Standard deviation of per-fragment signal ratio deviations.
        /// True peptides: low (all fragments deviate similarly, near 0). Interference: high.
        /// </summary>
        public float StdSignalRatioDeviation { get; set; }

        #endregion

        #region Category 4: Peak Shape

        /// <summary>Cosine at the peak-detected apex (may differ from full-window apex).</summary>
        public float PeakApexScore { get; set; }

        /// <summary>Temporal cosine restricted to detected peak boundaries.</summary>
        public float PeakTemporalScore { get; set; }

        /// <summary>
        /// Ratio of average boundary signal to apex signal.
        /// True peaks: low (signal drops at boundaries). Interference: high (doesn't drop).
        /// LOWER is better for targets.
        /// </summary>
        public float BoundarySignalRatio { get; set; }

        /// <summary>
        /// Ratio of apex signal to mean signal across the peak/window.
        /// True peaks: > 1 (apex stands above average). Flat noise: ≈ 1.
        /// HIGHER is better for targets.
        /// </summary>
        public float ApexToMeanRatio { get; set; }

        /// <summary>
        /// Number of candidate peak groups detected. More peaks in the window may
        /// indicate multiple co-eluting species (ambiguity).
        /// </summary>
        public int CandidateCount { get; set; }

        /// <summary>Detected peak group from DiaPeakGroupDetector. Null if no peak detected.</summary>
        public PeakGroup? DetectedPeakGroup { get; set; }

        /// <summary>Peak width in minutes (from detected peak boundaries).</summary>
        public float PeakWidth { get; set; }

        #endregion

        #region Category 5: Retention Time Quality

        /// <summary>Observed apex RT from extraction (minutes).</summary>
        public float ObservedApexRt { get; set; }

        /// <summary>Time point index of the observed apex.</summary>
        public int ApexTimeIndex { get; set; }

        /// <summary>Scoring strategy actually used for this result.</summary>
        public ScoringStrategy ScoringStrategyUsed { get; set; }

        #endregion

        #region Fragment Evidence

        /// <summary>Number of library fragment ions that yielded at least one XIC data point</summary>
        public int FragmentsDetected { get; set; }

        /// <summary>Total number of library fragment ions queried</summary>
        public int FragmentsQueried { get; }

        /// <summary>
        /// Sum of extracted intensities per fragment (length = FragmentsQueried).
        /// Parallel to the library's fragment ion order.
        /// </summary>
        public float[] ExtractedIntensities { get; }

        /// <summary>
        /// Number of XIC data points per fragment (length = FragmentsQueried).
        /// </summary>
        public int[] XicPointCounts { get; }

        /// <summary>Number of time points used for temporal scoring.</summary>
        public int TimePointsUsed { get; set; }

        #endregion

        #region Retention Time Context

        /// <summary>Library/predicted retention time (minutes). Null if library had no RT.</summary>
        public double? LibraryRetentionTime { get; }

        /// <summary>Lower bound of the RT extraction window (minutes)</summary>
        public float RtWindowStart { get; }

        /// <summary>Upper bound of the RT extraction window (minutes)</summary>
        public float RtWindowEnd { get; }

        #endregion

        #region Classifier & FDR

        /// <summary>
        /// Combined score from the linear discriminant / logistic regression classifier.
        /// Set by DiaLinearDiscriminant after feature extraction. Used for ranking and FDR.
        /// </summary>
        public float ClassifierScore { get; set; }

        /// <summary>
        /// FDR information computed by DiaFdrEngine (q-value, cumulative counts, etc.).
        /// </summary>
        public DiaFdrInfo FdrInfo { get; set; }

        #endregion

        public DiaSearchResult(
            string sequence,
            int chargeState,
            double precursorMz,
            int windowId,
            bool isDecoy,
            int fragmentsQueried,
            double? libraryRetentionTime,
            float rtWindowStart,
            float rtWindowEnd)
        {
            Sequence = sequence ?? throw new ArgumentNullException(nameof(sequence));
            ChargeState = chargeState;
            PrecursorMz = precursorMz;
            WindowId = windowId;
            IsDecoy = isDecoy;
            FragmentsQueried = fragmentsQueried;
            ExtractedIntensities = new float[fragmentsQueried];
            XicPointCounts = new int[fragmentsQueried];
            LibraryRetentionTime = libraryRetentionTime;
            RtWindowStart = rtWindowStart;
            RtWindowEnd = rtWindowEnd;

            // Initialize all scores to NaN
            DotProductScore = float.NaN;
            SpectralAngleScore = float.NaN;
            RawCosine = float.NaN;
            ApexDotProductScore = float.NaN;
            TemporalCosineScore = float.NaN;
            MeanFragmentCorrelation = float.NaN;
            MinFragmentCorrelation = float.NaN;
            PeakMeanFragCorrelation = float.NaN;
            PeakMinFragCorrelation = float.NaN;
            BestFragCorrelationSum = float.NaN;
            MedianFragRefCorr = float.NaN;
            MinFragRefCorr = float.NaN;
            StdFragRefCorr = float.NaN;
            BestFragWeightedCosine = float.NaN;
            BestFragIndex = -1;
            SmoothedMeanFragCorr = float.NaN;
            SmoothedMinFragCorr = float.NaN;
            LogTotalIntensity = 0f;
            Log2SignalToNoise = float.NaN;
            MeanSignalRatioDeviation = float.NaN;
            MaxSignalRatioDeviation = float.NaN;
            StdSignalRatioDeviation = float.NaN;
            PeakApexScore = float.NaN;
            PeakTemporalScore = float.NaN;
            BoundarySignalRatio = float.NaN;
            ApexToMeanRatio = float.NaN;
            CandidateCount = 0;
            PeakWidth = 0f;
            ClassifierScore = float.NaN;
            FdrInfo = null;
        }

        /// <summary>Whether this result meets the minimum fragment detection threshold.</summary>
        public bool MeetsMinFragments(int minRequired) => FragmentsDetected >= minRequired;

        /// <summary>Fraction of queried fragments that were detected (0–1).</summary>
        public float FragmentDetectionRate =>
            FragmentsQueried > 0 ? (float)FragmentsDetected / FragmentsQueried : 0f;

        /// <summary>
        /// RT deviation from predicted in minutes (signed).
        /// Positive = observed later than predicted.
        /// </summary>
        public float RtDeviationMinutes =>
            LibraryRetentionTime.HasValue ? ObservedApexRt - (float)LibraryRetentionTime.Value : 0f;

        /// <summary>
        /// Squared RT deviation (always positive). Used as a penalty feature.
        /// </summary>
        public float RtDeviationSquared => RtDeviationMinutes * RtDeviationMinutes;

        /// <summary>
        /// Intensity coefficient of variation across fragments at apex.
        /// Low CV = consistent intensities. High CV = variable (may indicate interference).
        /// </summary>
        public float IntensityCV
        {
            get
            {
                if (FragmentsDetected < 2) return float.NaN;
                float sum = 0f, sum2 = 0f;
                int n = 0;
                for (int i = 0; i < FragmentsQueried; i++)
                {
                    if (ExtractedIntensities[i] > 0f)
                    {
                        sum += ExtractedIntensities[i];
                        sum2 += ExtractedIntensities[i] * ExtractedIntensities[i];
                        n++;
                    }
                }
                if (n < 2) return float.NaN;
                float mean = sum / n;
                float variance = (sum2 / n) - (mean * mean);
                return mean > 0f ? MathF.Sqrt(Math.Max(0f, variance)) / mean : float.NaN;
            }
        }

        /// <summary>
        /// Builds the feature vector for the linear discriminant classifier.
        /// 
        /// Phase 13 recommended feature set (~25 features):
        /// - Drops PeakApexScore (inflates decoys, Phase 12 finding)
        /// - Drops PeakSymmetry (zero separation, Phase 12 finding)
        /// - Replaces MeanFragCorr with PeakMeanFragCorr (gap 0.47 → 0.69)
        /// - Adds best-fragment, signal ratio, smoothed correlation, S/N, peak shape features
        /// 
        /// Returns NaN-safe values (NaN → 0 for the classifier).
        /// </summary>
        public float[] GetFeatureVector()
        {
            return new float[]
            {
                // Category 1: Spectral Similarity (4 features)
                Safe(ApexDotProductScore),          // 0: Cosine at full-window apex
                Safe(TemporalCosineScore),          // 1: Intensity-weighted cosine across window
                Safe(SpectralAngleScore),           // 2: Spectral angle
                Safe(BestFragWeightedCosine),       // 3: Reference-weighted cosine at apex

                // Category 2: Fragment Co-elution (7 features)
                Safe(PeakMeanFragCorrelation),      // 4: Mean corr within peak (replaces full-window)
                Safe(PeakMinFragCorrelation),       // 5: Min corr within peak
                Safe(BestFragCorrelationSum),        // 6: Best-fragment total correlation sum
                Safe(MedianFragRefCorr),            // 7: Median fragment-vs-reference correlation
                Safe(MinFragRefCorr),               // 8: Min fragment-vs-reference correlation
                Safe(StdFragRefCorr),               // 9: Std of fragment-vs-reference correlations
                Safe(SmoothedMeanFragCorr),         // 10: Smoothed mean correlation

                // Category 3: Signal/Noise and Intensity (5 features)
                Safe(LogTotalIntensity),            // 11: Log2 total intensity
                Safe(Log2SignalToNoise),            // 12: Log2 signal-to-noise
                Safe(MeanSignalRatioDeviation),     // 13: Mean signal ratio deviation (LOWER=better)
                Safe(MaxSignalRatioDeviation),      // 14: Max signal ratio deviation (LOWER=better)
                Safe(StdSignalRatioDeviation),      // 15: Std signal ratio deviation

                // Category 4: Peak Shape (4 features)
                Safe(BoundarySignalRatio),          // 16: Boundary/apex ratio (LOWER=better)
                Safe(ApexToMeanRatio),              // 17: Apex/mean ratio (HIGHER=better)
                (float)CandidateCount,              // 18: Number of peak candidates
                Safe(PeakWidth),                    // 19: Peak width in minutes

                // Category 5: RT Quality (2 features)
                Safe(RtDeviationMinutes),           // 20: RT deviation (signed)
                Safe(RtDeviationSquared),           // 21: RT deviation squared

                // Category 6: Fragment Evidence (1 feature)
                FragmentDetectionRate,              // 22: Fraction of fragments detected
            };
        }

        /// <summary>
        /// Number of features in the feature vector.
        /// </summary>
        public static int FeatureCount => 23;

        /// <summary>
        /// Feature names for reporting and diagnostics.
        /// </summary>
        public static string[] FeatureNames => new[]
        {
            "ApexScore", "TemporalScore", "SpectralAngle", "BestFragWeightedCosine",
            "PeakMeanFragCorr", "PeakMinFragCorr", "BestFragCorrSum", "MedianFragRefCorr",
            "MinFragRefCorr", "StdFragRefCorr", "SmoothedMeanFragCorr",
            "LogTotalIntensity", "Log2SNR", "MeanSigRatioDev", "MaxSigRatioDev", "StdSigRatioDev",
            "BoundarySignalRatio", "ApexToMeanRatio", "CandidateCount", "PeakWidth",
            "RtDevMinutes", "RtDevSquared",
            "FragDetRate"
        };

        /// <summary>Replace NaN with 0 for classifier input.</summary>
        private static float Safe(float v) => float.IsNaN(v) ? 0f : v;

        public override string ToString()
        {
            return $"{Sequence}/{ChargeState} Window={WindowId} " +
                   $"DotProduct={DotProductScore:F4} SpectralAngle={SpectralAngleScore:F4} " +
                   $"BestFragCorrSum={BestFragCorrelationSum:F3} MeanSigRatioDev={MeanSignalRatioDeviation:F3} " +
                   $"Fragments={FragmentsDetected}/{FragmentsQueried}" +
                   (IsDecoy ? " [DECOY]" : "");
        }
    }
}