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
    /// </summary>
    public struct DiaFeatureVector
    {
        // ── Primary scores (from Phase 9 temporal scoring) ──────────────────

        /// <summary>Cosine at the consensus chromatographic apex [0,1]</summary>
        public float ApexScore;

        /// <summary>Average cosine across time points [0,1]</summary>
        public float TemporalScore;

        /// <summary>Raw cosine before any nonlinear transform [0,1]</summary>
        public float RawCosine;

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

        /// <summary>RT window half-width used for extraction (minutes).</summary>
        public float RtWindowHalfWidth;

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
        /// Updated from 12 to 14 with coelution features.
        /// </summary>
        public const int ClassifierFeatureCount = 14;

        /// <summary>
        /// Writes the classifier features into a float span for linear algebra operations.
        /// Order must be consistent with weight vectors.
        /// </summary>
        public readonly void WriteTo(Span<float> features)
        {
            if (features.Length < ClassifierFeatureCount)
                throw new ArgumentException($"Span must have at least {ClassifierFeatureCount} elements");

            features[0] = ApexScore;
            features[1] = TemporalScore;
            features[2] = RawCosine;
            features[3] = SpectralAngle;
            features[4] = MeanFragmentCorrelation;
            features[5] = MinFragmentCorrelation;
            features[6] = FragmentDetectionRate;
            features[7] = LogTotalIntensity;
            features[8] = IntensityCV;
            features[9] = MedianXicDepth;
            features[10] = XicDepthCV;
            features[11] = (float)TimePointsUsed;
            features[12] = RtDeviationMinutes;
            features[13] = RtWindowHalfWidth;
        }

        /// <summary>Feature names in the same order as WriteTo, for reporting.</summary>
        public static readonly string[] FeatureNames = new[]
        {
            "ApexScore",
            "TemporalScore",
            "RawCosine",
            "SpectralAngle",
            "MeanFragCorr",
            "MinFragCorr",
            "FragmentDetectionRate",
            "LogTotalIntensity",
            "IntensityCV",
            "MedianXicDepth",
            "XicDepthCV",
            "TimePointsUsed",
            "RtDeviationMinutes",
            "RtWindowHalfWidth",
        };
    }

    /// <summary>
    /// Computes feature vectors from DiaSearchResult objects.
    /// Thread-safe: uses only local state + ArrayPool.
    /// </summary>
    public static class DiaFeatureExtractor
    {
        /// <summary>
        /// Computes a feature vector from a scored DiaSearchResult.
        /// Now uses ObservedApexRt for real RT deviation and
        /// MeanFragmentCorrelation/MinFragmentCorrelation from the result.
        /// </summary>
        public static DiaFeatureVector ComputeFeatures(
            DiaSearchResult result,
            int precursorIndex)
        {
            var fv = new DiaFeatureVector();

            // ── Primary scores ──────────────────────────────────────────
            fv.ApexScore = SafeScore(result.ApexDotProductScore);
            fv.TemporalScore = SafeScore(result.TemporalCosineScore);
            fv.RawCosine = SafeScore(result.RawCosine);
            fv.SpectralAngle = SafeScore(result.SpectralAngleScore);

            // ── Coelution features (new in Phase 10 revision) ───────────
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

            // ── Retention time features (now using ObservedApexRt) ──────
            fv.RtWindowHalfWidth = (result.RtWindowEnd - result.RtWindowStart) / 2f;

            if (result.LibraryRetentionTime.HasValue && !float.IsNaN(result.ObservedApexRt))
            {
                fv.RtDeviationMinutes = MathF.Abs(
                    result.ObservedApexRt - (float)result.LibraryRetentionTime.Value);
            }
            else
            {
                fv.RtDeviationMinutes = fv.RtWindowHalfWidth; // fallback
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