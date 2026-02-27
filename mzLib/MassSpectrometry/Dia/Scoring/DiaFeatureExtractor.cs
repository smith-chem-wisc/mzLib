// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Buffers;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Feature vector for a single DIA precursor identification.
    /// Contains all features that can feed into a discriminant classifier.
    /// 
    /// Designed as a struct with fixed-size fields for cache-friendly batch processing.
    /// All features are normalized to roughly comparable scales where possible.
    /// 
    /// Phase 10.5 changes:
    ///   - REMOVED RawCosine (r=1.000 with TemporalScore — identical, caused LDA degeneracy)
    ///   - REMOVED RtWindowHalfWidth (constant across all precursors — zero variance, zero info)
    ///   - ADDED RtDeviationSquared: (ΔRT)² quadratic penalty for RT outliers (DIA-NN approach)
    /// Phase 12 changes:
    ///   - ADDED PeakApexScore: cosine at peak-detected apex (more robust than max-signal apex)
    ///   - ADDED PeakMeanFragCorr: fragment correlation within peak boundaries only
    ///   - ADDED PeakWidth: chromatographic peak width in minutes
    ///   - ADDED PeakSymmetry: peak shape symmetry ratio
    /// Feature count: 14 → 12 → 13 → 17 (after Phase 12 peak group features)
    /// </summary>
    public struct DiaFeatureVector
    {
        // ── Primary scores (from Phase 9 temporal scoring) ──────────────────

        /// <summary>Cosine at the consensus chromatographic apex [0,1]</summary>
        public float ApexScore;

        /// <summary>Average cosine across time points [0,1]</summary>
        public float TemporalScore;

        // RawCosine REMOVED in Phase 10.5a (r=1.000 with TemporalScore — identical)

        /// <summary>Spectral angle score [0,1]</summary>
        public float SpectralAngle;

        // ── Coelution features (Phase 10 — highest-leverage additions) ──────

        /// <summary>
        /// Mean pairwise Pearson correlation across detected fragment XICs.
        /// True coeluting fragments: ~0.9+. Interfered: much lower.
        /// This is the single most discriminative non-score feature per DIA-NN.
        /// </summary>
        public float MeanFragmentCorrelation;

        /// <summary>
        /// Minimum pairwise Pearson correlation — the weakest link.
        /// Identifies the most interfered fragment.
        /// </summary>
        public float MinFragmentCorrelation;

        // ── Fragment evidence features ──────────────────────────────────────

        /// <summary>Fraction of queried fragments that had XIC data [0,1]</summary>
        public float FragmentDetectionRate;

        // ── Intensity features ──────────────────────────────────────────────

        /// <summary>Log10 of total extracted intensity across all fragments.</summary>
        public float LogTotalIntensity;

        /// <summary>CV of fragment intensities across detected fragments.</summary>
        public float IntensityCV;

        // ── XIC shape features ──────────────────────────────────────────────

        /// <summary>Median number of XIC data points per detected fragment.</summary>
        public float MedianXicDepth;

        /// <summary>CV of XIC point counts across detected fragments.</summary>
        public float XicDepthCV;

        // ── Temporal evidence ───────────────────────────────────────────────

        /// <summary>Number of RT time points that contributed to the temporal score.</summary>
        public int TimePointsUsed;

        // ── Retention time features ─────────────────────────────────────────

        /// <summary>|Observed apex RT - Library RT| in minutes. Lower = better.</summary>
        public float RtDeviationMinutes;

        /// <summary>
        /// (ΔRT)² — squared RT deviation in minutes².
        /// Quadratic penalty that sharpens discrimination for large RT outliers.
        /// DIA-NN explicitly uses this form. Complements the linear RtDeviationMinutes.
        /// Added in Phase 10.5b.
        /// </summary>
        public float RtDeviationSquared;

        // RtWindowHalfWidth REMOVED in Phase 10.5a (constant, zero information)

        // ── Peak group features (Phase 12) ──────────────────────────────

        /// <summary>
        /// Cosine similarity at the peak-group-detected apex [0,1].
        /// Unlike ApexScore (max-signal scan, which can be an interference spike),
        /// this uses the apex from proper peak boundary detection.
        /// Falls back to ApexScore if no peak was detected.
        /// </summary>
        public float PeakApexScore;

        /// <summary>
        /// Mean pairwise fragment correlation computed ONLY within the detected
        /// peak boundaries. By restricting to the actual elution peak, interference
        /// from flanking co-eluters is excluded, dramatically improving discrimination.
        /// Falls back to MeanFragmentCorrelation if no peak was detected.
        /// </summary>
        public float PeakMeanFragCorr;

        /// <summary>
        /// Chromatographic peak width in minutes. Narrower peaks indicate sharper,
        /// more confident detections. Very wide "peaks" may indicate interference
        /// or low-abundance signal that doesn't form a discrete peak.
        /// 0 if no peak was detected.
        /// </summary>
        public float PeakWidth;

        /// <summary>
        /// Peak shape symmetry ratio: (apex - left) / (right - left).
        /// 0.5 = perfectly symmetric. Real chromatographic peaks are typically 0.3–0.7.
        /// Values near 0 or 1 indicate truncated peaks (edge of RT window).
        /// 0.5 if no peak was detected.
        /// </summary>
        public float PeakSymmetry;

        // ── Metadata (not classifier features) ─────────────────────────────

        /// <summary>Whether this is a decoy precursor</summary>
        public bool IsDecoy;

        /// <summary>Index into the original result list for traceability</summary>
        public int PrecursorIndex;

        /// <summary>Precursor charge state</summary>
        public int ChargeState;

        /// <summary>Precursor m/z</summary>
        public float PrecursorMz;

        /// <summary>Number of fragments detected</summary>
        public int FragmentsDetected;

        /// <summary>Number of fragments queried</summary>
        public int FragmentsQueried;

        /// <summary>
        /// Number of features used by the classifier.
        /// Phase 12: 13 → 17 (added PeakApexScore, PeakMeanFragCorr, PeakWidth, PeakSymmetry).
        /// </summary>
        public const int ClassifierFeatureCount = 17;

        /// <summary>
        /// Writes the classifier features into a float span for linear algebra operations.
        /// Order must be consistent with weight vectors.
        /// 
        /// Phase 12 order (17 features):
        ///   [0] ApexScore, [1] TemporalScore, [2] SpectralAngle,
        ///   [3] MeanFragCorr, [4] MinFragCorr, [5] FragDetRate,
        ///   [6] LogTotalIntensity, [7] IntensityCV, [8] MedianXicDepth,
        ///   [9] XicDepthCV, [10] TimePointsUsed, [11] RtDeviationMinutes,
        ///   [12] RtDeviationSquared, [13] PeakApexScore, [14] PeakMeanFragCorr,
        ///   [15] PeakWidth, [16] PeakSymmetry
        /// </summary>
        public readonly void WriteTo(Span<float> features)
        {
            if (features.Length < ClassifierFeatureCount)
                throw new ArgumentException($"Span must have at least {ClassifierFeatureCount} elements");

            features[0] = ApexScore;
            features[1] = TemporalScore;
            features[2] = SpectralAngle;
            features[3] = MeanFragmentCorrelation;
            features[4] = MinFragmentCorrelation;
            features[5] = FragmentDetectionRate;
            features[6] = LogTotalIntensity;
            features[7] = IntensityCV;
            features[8] = MedianXicDepth;
            features[9] = XicDepthCV;
            features[10] = (float)TimePointsUsed;
            features[11] = RtDeviationMinutes;
            features[12] = RtDeviationSquared;
            features[13] = PeakApexScore;
            features[14] = PeakMeanFragCorr;
            features[15] = PeakWidth;
            features[16] = PeakSymmetry;
        }

        /// <summary>Feature names in the same order as WriteTo, for reporting.</summary>
        public static readonly string[] FeatureNames = new[]
        {
            "ApexScore",
            "TemporalScore",
            "SpectralAngle",
            "MeanFragCorr",
            "MinFragCorr",
            "FragDetRate",
            "LogTotalIntensity",
            "IntensityCV",
            "MedianXicDepth",
            "XicDepthCV",
            "TimePointsUsed",
            "RtDeviationMinutes",
            "RtDeviationSquared",
            "PeakApexScore",
            "PeakMeanFragCorr",
            "PeakWidth",
            "PeakSymmetry",
        };
    }

    /// <summary>
    /// Computes feature vectors from DiaSearchResult objects.
    /// Thread-safe: uses only local state + ArrayPool.
    /// 
    /// Phase 10.5 changes:
    ///   - Removed RawCosine population (was identical to TemporalScore)
    ///   - Removed RtWindowHalfWidth computation (constant, zero information)
    ///   - Added RtDeviationSquared = (ΔRT)² computation
    /// </summary>
    public static class DiaFeatureExtractor
    {
        /// <summary>
        /// Computes a feature vector from a scored DiaSearchResult.
        /// </summary>
        public static DiaFeatureVector ComputeFeatures(
            DiaSearchResult result,
            int precursorIndex)
        {
            var fv = new DiaFeatureVector();

            // ── Primary scores ──────────────────────────────────────────
            fv.ApexScore = SafeScore(result.ApexDotProductScore);
            fv.TemporalScore = SafeScore(result.TemporalCosineScore);
            // RawCosine REMOVED — was identical to TemporalScore (r=1.000)
            fv.SpectralAngle = SafeScore(result.SpectralAngleScore);

            // ── Coelution features ──────────────────────────────────────
            fv.MeanFragmentCorrelation = SafeScore(result.MeanFragmentCorrelation);
            fv.MinFragmentCorrelation = float.IsNaN(result.MinFragmentCorrelation)
                ? -1f  // worst case sentinel for NaN
                : result.MinFragmentCorrelation;

            // ── Fragment evidence ───────────────────────────────────────
            fv.FragmentDetectionRate = result.FragmentDetectionRate;
            fv.FragmentsDetected = result.FragmentsDetected;
            fv.FragmentsQueried = result.FragmentsQueried;

            // ── Temporal evidence ───────────────────────────────────────
            fv.TimePointsUsed = result.TimePointsUsed;

            // ── Intensity features ─────────────────────────────────────
            float totalIntensity = 0f;
            int nDetected = 0;

            for (int f = 0; f < result.FragmentsQueried; f++)
            {
                totalIntensity += result.ExtractedIntensities[f];
                if (result.XicPointCounts[f] > 0)
                    nDetected++;
            }

            fv.LogTotalIntensity = totalIntensity > 0 ? MathF.Log10(totalIntensity) : 0f;

            if (nDetected >= 2)
            {
                float mean = 0f;
                for (int f = 0; f < result.FragmentsQueried; f++)
                    if (result.XicPointCounts[f] > 0)
                        mean += result.ExtractedIntensities[f];
                mean /= nDetected;

                float variance = 0f;
                for (int f = 0; f < result.FragmentsQueried; f++)
                {
                    if (result.XicPointCounts[f] > 0)
                    {
                        float diff = result.ExtractedIntensities[f] - mean;
                        variance += diff * diff;
                    }
                }
                variance /= nDetected;
                fv.IntensityCV = mean > 0 ? MathF.Sqrt(variance) / mean : 0f;
            }
            else
            {
                fv.IntensityCV = 1f;
            }

            // ── XIC shape features ─────────────────────────────────────
            if (nDetected >= 1)
            {
                int[] counts = ArrayPool<int>.Shared.Rent(nDetected);
                try
                {
                    int idx = 0;
                    for (int f = 0; f < result.FragmentsQueried; f++)
                        if (result.XicPointCounts[f] > 0)
                            counts[idx++] = result.XicPointCounts[f];

                    Array.Sort(counts, 0, nDetected);
                    fv.MedianXicDepth = nDetected % 2 == 1
                        ? counts[nDetected / 2]
                        : (counts[nDetected / 2 - 1] + counts[nDetected / 2]) / 2f;

                    if (nDetected >= 2)
                    {
                        float mean = 0f;
                        for (int i = 0; i < nDetected; i++) mean += counts[i];
                        mean /= nDetected;
                        float var2 = 0f;
                        for (int i = 0; i < nDetected; i++)
                        {
                            float diff = counts[i] - mean;
                            var2 += diff * diff;
                        }
                        var2 /= nDetected;
                        fv.XicDepthCV = mean > 0 ? MathF.Sqrt(var2) / mean : 0f;
                    }
                    else
                        fv.XicDepthCV = 1f;
                }
                finally
                {
                    ArrayPool<int>.Shared.Return(counts);
                }
            }
            else
            {
                fv.MedianXicDepth = 0f;
                fv.XicDepthCV = 1f;
            }

            // ── Retention time features ─────────────────────────────────
            // Maximum RT deviation cap to prevent pathological scores.
            // Without this, precursors with no library RT get rtWindowHalfWidth
            // (which can be ~22 min for full-run fallback), causing
            // RtDeviationSquared ≈ 484 and catastrophic classifier scores.
            const float MaxRtDeviationMinutes = 5.0f;

            if (result.LibraryRetentionTime.HasValue && !float.IsNaN(result.ObservedApexRt))
            {
                float deltaRt = MathF.Abs(
                    result.ObservedApexRt - (float)result.LibraryRetentionTime.Value);
                deltaRt = MathF.Min(deltaRt, MaxRtDeviationMinutes);
                fv.RtDeviationMinutes = deltaRt;
                fv.RtDeviationSquared = deltaRt * deltaRt;
            }
            else
            {
                // No library RT or no observed apex — use maximum penalty
                fv.RtDeviationMinutes = MaxRtDeviationMinutes;
                fv.RtDeviationSquared = MaxRtDeviationMinutes * MaxRtDeviationMinutes;
            }

            // ── Peak group features (Phase 12) ──────────────────────────
            if (result.DetectedPeakGroup.HasValue && result.DetectedPeakGroup.Value.IsValid)
            {
                var pg = result.DetectedPeakGroup.Value;

                fv.PeakApexScore = SafeScore(result.PeakApexScore);
                fv.PeakMeanFragCorr = SafeScore(result.PeakMeanFragCorrelation);
                fv.PeakWidth = pg.PeakWidthMinutes;
                fv.PeakSymmetry = pg.SymmetryRatio;
            }
            else
            {
                // Fall back to full-window scores when no peak was detected
                fv.PeakApexScore = SafeScore(result.ApexDotProductScore);
                fv.PeakMeanFragCorr = SafeScore(result.MeanFragmentCorrelation);
                fv.PeakWidth = 0f;
                fv.PeakSymmetry = 0.5f;
            }

            // ── Metadata ────────────────────────────────────────────────
            fv.ChargeState = result.ChargeState;
            fv.PrecursorMz = (float)result.PrecursorMz;
            fv.IsDecoy = result.IsDecoy;
            fv.PrecursorIndex = precursorIndex;

            return fv;
        }

        /// <summary>Replaces NaN scores with 0 for classifier safety.</summary>
        private static float SafeScore(float score) =>
            float.IsNaN(score) ? 0f : score;
    }
}