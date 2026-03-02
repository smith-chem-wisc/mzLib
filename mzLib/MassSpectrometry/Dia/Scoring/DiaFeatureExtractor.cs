// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Buffers;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Feature vector for a single DIA precursor identification.
    /// 
    /// Phase 13: 17 → 28 features. Adds best-fragment reference curve,
    /// signal ratio deviation, smoothed correlations, S/N, and peak shape.
    /// 
    /// Drops PeakSymmetry (separation 0.01 — zero discriminative value).
    /// </summary>
    public struct DiaFeatureVector
    {
        // ── Primary scores ──────────────────────────────────────────
        public float ApexScore;
        public float TemporalScore;
        public float SpectralAngle;

        // ── Coelution features (Phase 10) ───────────────────────────
        public float MeanFragmentCorrelation;
        public float MinFragmentCorrelation;

        // ── Fragment evidence ───────────────────────────────────────
        public float FragmentDetectionRate;

        // ── Intensity features ─────────────────────────────────────
        public float LogTotalIntensity;
        public float IntensityCV;

        // ── XIC shape features ─────────────────────────────────────
        public float MedianXicDepth;
        public float XicDepthCV;

        // ── Temporal evidence ──────────────────────────────────────
        public int TimePointsUsed;

        // ── Retention time ─────────────────────────────────────────
        public float RtDeviationMinutes;
        public float RtDeviationSquared;

        // ── Peak group features (Phase 12) ─────────────────────────
        public float PeakApexScore;
        public float PeakMeanFragCorr;
        public float PeakWidth;
        // PeakSymmetry DROPPED: separation = 0.01, zero discriminative value

        // ── Phase 13: Best-fragment reference curve ─────────────────
        public float BestFragCorrelationSum;
        public float MedianFragRefCorr;
        public float MinFragRefCorr;
        public float StdFragRefCorr;
        public float BestFragWeightedCosine;

        // ── Phase 13: Smoothed correlations ─────────────────────────
        public float SmoothedMeanFragCorr;

        // ── Phase 13: Signal ratio deviation ────────────────────────
        public float MeanSignalRatioDeviation;
        public float MaxSignalRatioDeviation;
        public float StdSignalRatioDeviation;

        // ── Phase 13: Signal-to-noise ───────────────────────────────
        public float Log2SignalToNoise;

        // ── Phase 13: Peak shape ────────────────────────────────────
        public float BoundarySignalRatio;
        public float ApexToMeanRatio;

        // ── Metadata (not classifier features) ─────────────────────
        public bool IsDecoy;
        public int PrecursorIndex;
        public int ChargeState;
        public float PrecursorMz;
        public int FragmentsDetected;
        public int FragmentsQueried;

        /// <summary>
        /// Number of features used by the classifier.
        /// Phase 13: 17 → 28 (dropped PeakSymmetry, added 12 new features).
        /// </summary>
        public const int ClassifierFeatureCount = 28;

        /// <summary>
        /// Writes classifier features into a float span.
        /// Order must be consistent with weight vectors.
        /// 
        /// Phase 13 order (28 features):
        ///   [0-2]   Primary scores: ApexScore, TemporalScore, SpectralAngle
        ///   [3-4]   Coelution: MeanFragCorr, MinFragCorr
        ///   [5]     Fragment evidence: FragDetRate
        ///   [6-7]   Intensity: LogTotalIntensity, IntensityCV
        ///   [8-9]   XIC shape: MedianXicDepth, XicDepthCV
        ///   [10]    Temporal: TimePointsUsed
        ///   [11-12] RT: RtDeviationMinutes, RtDeviationSquared
        ///   [13-15] Peak group: PeakApexScore, PeakMeanFragCorr, PeakWidth
        ///   [16-20] Best-fragment: BestFragCorrSum, MedianRefCorr, MinRefCorr, StdRefCorr, WeightedCosine
        ///   [21]    Smoothed: SmoothedMeanFragCorr
        ///   [22-24] Signal ratio: MeanSigRatioDev, MaxSigRatioDev, StdSigRatioDev
        ///   [25]    S/N: Log2SignalToNoise
        ///   [26-27] Peak shape: BoundarySignalRatio, ApexToMeanRatio
        /// </summary>
        public readonly void WriteTo(Span<float> features)
        {
            if (features.Length < ClassifierFeatureCount)
                throw new ArgumentException($"Span must have at least {ClassifierFeatureCount} elements");

            // Existing features [0-15] — same order as Phase 12 minus PeakSymmetry
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

            // Phase 13 features [16-27]
            features[16] = BestFragCorrelationSum;
            features[17] = MedianFragRefCorr;
            features[18] = MinFragRefCorr;
            features[19] = StdFragRefCorr;
            features[20] = BestFragWeightedCosine;
            features[21] = SmoothedMeanFragCorr;
            features[22] = MeanSignalRatioDeviation;
            features[23] = MaxSignalRatioDeviation;
            features[24] = StdSignalRatioDeviation;
            features[25] = Log2SignalToNoise;
            features[26] = BoundarySignalRatio;
            features[27] = ApexToMeanRatio;
        }

        public static readonly string[] FeatureNames = new[]
        {
            "ApexScore", "TemporalScore", "SpectralAngle",
            "MeanFragCorr", "MinFragCorr", "FragDetRate",
            "LogTotalIntensity", "IntensityCV",
            "MedianXicDepth", "XicDepthCV",
            "TimePointsUsed",
            "RtDeviationMinutes", "RtDeviationSquared",
            "PeakApexScore", "PeakMeanFragCorr", "PeakWidth",
            "BestFragCorrSum", "MedianFragRefCorr", "MinFragRefCorr",
            "StdFragRefCorr", "BestFragWeightedCosine",
            "SmoothedMeanFragCorr",
            "MeanSigRatioDev", "MaxSigRatioDev", "StdSigRatioDev",
            "Log2SNR",
            "BoundarySignalRatio", "ApexToMeanRatio",
        };
    }

    /// <summary>
    /// Computes feature vectors from DiaSearchResult objects.
    /// Thread-safe: uses only local state + ArrayPool.
    /// 
    /// Phase 13: populates 28-feature vector including best-fragment,
    /// signal ratio, smoothed correlation, S/N, and peak shape features.
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

            // ── Primary scores ──────────────────────────────────────
            fv.ApexScore = SafeScore(result.ApexDotProductScore);
            fv.TemporalScore = SafeScore(result.TemporalCosineScore);
            fv.SpectralAngle = SafeScore(result.SpectralAngleScore);

            // ── Coelution features ──────────────────────────────────
            fv.MeanFragmentCorrelation = SafeScore(result.MeanFragmentCorrelation);
            fv.MinFragmentCorrelation = float.IsNaN(result.MinFragmentCorrelation)
                ? -1f
                : result.MinFragmentCorrelation;

            // ── Fragment evidence ───────────────────────────────────
            fv.FragmentDetectionRate = result.FragmentDetectionRate;
            fv.FragmentsDetected = result.FragmentsDetected;
            fv.FragmentsQueried = result.FragmentsQueried;

            // ── Temporal evidence ───────────────────────────────────
            fv.TimePointsUsed = result.TimePointsUsed;

            // ── Intensity features ─────────────────────────────────
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

            // ── XIC shape features ─────────────────────────────────
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

            // ── Retention time features ─────────────────────────────
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
                fv.RtDeviationMinutes = MaxRtDeviationMinutes;
                fv.RtDeviationSquared = MaxRtDeviationMinutes * MaxRtDeviationMinutes;
            }

            // ── Peak group features (Phase 12) ──────────────────────
            if (result.DetectedPeakGroup.HasValue && result.DetectedPeakGroup.Value.IsValid)
            {
                var pg = result.DetectedPeakGroup.Value;

                fv.PeakApexScore = SafeScore(result.PeakApexScore);
                fv.PeakMeanFragCorr = SafeScore(result.PeakMeanFragCorrelation);
                fv.PeakWidth = pg.PeakWidthMinutes;
            }
            else
            {
                fv.PeakApexScore = SafeScore(result.ApexDotProductScore);
                fv.PeakMeanFragCorr = SafeScore(result.MeanFragmentCorrelation);
                fv.PeakWidth = 0f;
            }

            // ══════════════════════════════════════════════════════════
            //  Phase 13: New discriminative features
            //  Read directly from DiaSearchResult (populated by
            //  DiaFeatureCalculator during assembly)
            // ══════════════════════════════════════════════════════════

            // Best-fragment reference curve
            fv.BestFragCorrelationSum = SafeScore(result.BestFragCorrelationSum);
            fv.MedianFragRefCorr = SafeScore(result.MedianFragRefCorr);
            fv.MinFragRefCorr = float.IsNaN(result.MinFragRefCorr)
                ? -1f : result.MinFragRefCorr;
            fv.StdFragRefCorr = SafeScore(result.StdFragRefCorr);
            fv.BestFragWeightedCosine = SafeScore(result.BestFragWeightedCosine);

            // Smoothed correlations
            fv.SmoothedMeanFragCorr = SafeScore(result.SmoothedMeanFragCorr);

            // Signal ratio deviation (LOWER = better for targets)
            fv.MeanSignalRatioDeviation = SafeScore(result.MeanSignalRatioDeviation);
            fv.MaxSignalRatioDeviation = SafeScore(result.MaxSignalRatioDeviation);
            fv.StdSignalRatioDeviation = SafeScore(result.StdSignalRatioDeviation);

            // Signal-to-noise
            fv.Log2SignalToNoise = SafeScore(result.Log2SignalToNoise);

            // Peak shape
            fv.BoundarySignalRatio = SafeScore(result.BoundarySignalRatio);
            fv.ApexToMeanRatio = SafeScore(result.ApexToMeanRatio);

            // ── Metadata ────────────────────────────────────────────
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