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

        /// <summary>Cosine at the consensus chromatographic apex [0,1]. From DiaSearchResult.ApexDotProductScore.</summary>
        public float ApexScore;

        /// <summary>Average cosine across time points with sufficient fragments [0,1]. From DiaSearchResult.TemporalCosineScore.</summary>
        public float TemporalScore;

        /// <summary>Raw cosine before any nonlinear transform [0,1]. From DiaSearchResult.RawCosine.</summary>
        public float RawCosine;

        /// <summary>Spectral angle score [0,1]. From DiaSearchResult.SpectralAngleScore.</summary>
        public float SpectralAngle;

        // ── Fragment evidence features ──────────────────────────────────────

        /// <summary>Fraction of queried fragments that had XIC data [0,1]</summary>
        public float FragmentDetectionRate;

        // ── Intensity features ──────────────────────────────────────────────

        /// <summary>Log10 of total extracted intensity across all fragments. Higher = more signal.</summary>
        public float LogTotalIntensity;

        /// <summary>
        /// Coefficient of variation of fragment intensities (across detected fragments).
        /// Low CV = fragments have similar amounts of signal = consistent coelution.
        /// High CV = uneven fragment detection = possible interference.
        /// </summary>
        public float IntensityCV;

        // ── XIC shape features ──────────────────────────────────────────────

        /// <summary>
        /// Median number of XIC data points per detected fragment.
        /// More points = broader/better chromatographic coverage.
        /// </summary>
        public float MedianXicDepth;

        /// <summary>
        /// CV of XIC point counts across detected fragments.
        /// Low = fragments extracted consistently across similar number of scans.
        /// </summary>
        public float XicDepthCV;

        // ── Temporal evidence features ──────────────────────────────────────

        /// <summary>
        /// Number of RT time points that contributed to the temporal score.
        /// From DiaSearchResult.TimePointsUsed. Higher = more temporal evidence.
        /// </summary>
        public int TimePointsUsed;

        // ── Retention time features ─────────────────────────────────────────

        /// <summary>
        /// |Observed apex RT - Library RT| in minutes. Lower = better RT prediction.
        /// </summary>
        public float RtDeviationMinutes;

        /// <summary>
        /// RT window half-width used for extraction (minutes).
        /// </summary>
        public float RtWindowHalfWidth;

        // ── Metadata (not classifier features — for ground truth / tracing) ─

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
        /// </summary>
        public const int ClassifierFeatureCount = 12;

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
            features[4] = FragmentDetectionRate;
            features[5] = LogTotalIntensity;
            features[6] = IntensityCV;
            features[7] = MedianXicDepth;
            features[8] = XicDepthCV;
            features[9] = (float)TimePointsUsed;
            features[10] = RtDeviationMinutes;
            features[11] = RtWindowHalfWidth;
        }

        /// <summary>Feature names in the same order as WriteTo, for reporting.</summary>
        public static readonly string[] FeatureNames = new[]
        {
            "ApexScore",
            "TemporalScore",
            "RawCosine",
            "SpectralAngle",
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
    /// Computes feature vectors from DiaSearchResult objects and their XIC data.
    /// 
    /// This is the bridge between the extraction/scoring pipeline (Phases 1-9)
    /// and the classifier (Phase 10). It takes the existing per-precursor results
    /// and computes additional derived features useful for discriminating
    /// true from false identifications.
    /// 
    /// Thread-safe: uses only local state + ArrayPool.
    /// </summary>
    public static class DiaFeatureExtractor
    {
        /// <summary>
        /// Computes a feature vector from a scored DiaSearchResult.
        /// 
        /// Uses the Phase 9 hybrid scoring fields (ApexDotProductScore, TemporalCosineScore)
        /// plus derived features from the per-fragment intensity and XIC depth arrays.
        /// 
        /// When ExtractionResult and PrecursorQueryGroup are available, also computes
        /// the observed apex RT for RT deviation. Otherwise falls back to window half-width.
        /// </summary>
        public static DiaFeatureVector ComputeFeatures(
            DiaSearchResult result,
            int precursorIndex,
            ExtractionResult extractionResult = null,
            DiaLibraryQueryGenerator.PrecursorQueryGroup? queryGroup = null)
        {
            var fv = new DiaFeatureVector();

            // ── Primary scores (directly from Phase 9 temporal scorer) ──────
            fv.ApexScore = SafeScore(result.ApexDotProductScore);
            fv.TemporalScore = SafeScore(result.TemporalCosineScore);
            fv.RawCosine = SafeScore(result.RawCosine);
            fv.SpectralAngle = SafeScore(result.SpectralAngleScore);

            // ── Fragment evidence ───────────────────────────────────────────
            fv.FragmentDetectionRate = result.FragmentDetectionRate;
            fv.FragmentsDetected = result.FragmentsDetected;
            fv.FragmentsQueried = result.FragmentsQueried;

            // ── Temporal evidence ───────────────────────────────────────────
            fv.TimePointsUsed = result.TimePointsUsed;

            // ── Intensity features from per-fragment summed intensities ─────
            float totalIntensity = 0f;
            int nDetected = 0;

            for (int f = 0; f < result.FragmentsQueried; f++)
            {
                totalIntensity += result.ExtractedIntensities[f];
                if (result.XicPointCounts[f] > 0)
                    nDetected++;
            }

            fv.LogTotalIntensity = totalIntensity > 0 ? MathF.Log10(totalIntensity) : 0f;

            // Intensity CV across detected fragments
            if (nDetected >= 2)
            {
                float mean = 0f;
                for (int f = 0; f < result.FragmentsQueried; f++)
                {
                    if (result.XicPointCounts[f] > 0)
                        mean += result.ExtractedIntensities[f];
                }
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

            // ── XIC shape features from per-fragment point counts ───────────
            if (nDetected >= 1)
            {
                int[] counts = ArrayPool<int>.Shared.Rent(nDetected);
                try
                {
                    int idx = 0;
                    for (int f = 0; f < result.FragmentsQueried; f++)
                    {
                        if (result.XicPointCounts[f] > 0)
                            counts[idx++] = result.XicPointCounts[f];
                    }

                    // Median XIC depth
                    Array.Sort(counts, 0, nDetected);
                    fv.MedianXicDepth = nDetected % 2 == 1
                        ? counts[nDetected / 2]
                        : (counts[nDetected / 2 - 1] + counts[nDetected / 2]) / 2f;

                    // XIC depth CV
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
                    {
                        fv.XicDepthCV = 1f;
                    }
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

            // ── Retention time features ─────────────────────────────────────
            fv.RtWindowHalfWidth = (result.RtWindowEnd - result.RtWindowStart) / 2f;

            if (result.LibraryRetentionTime.HasValue &&
                extractionResult != null && queryGroup.HasValue &&
                result.ApexTimeIndex >= 0)
            {
                float observedApexRt = ComputeApexRt(
                    extractionResult, queryGroup.Value, result.ApexTimeIndex);
                if (!float.IsNaN(observedApexRt))
                {
                    fv.RtDeviationMinutes = MathF.Abs(
                        observedApexRt - (float)result.LibraryRetentionTime.Value);
                }
                else
                {
                    fv.RtDeviationMinutes = fv.RtWindowHalfWidth;
                }
            }
            else
            {
                fv.RtDeviationMinutes = fv.RtWindowHalfWidth;
            }

            // ── Metadata ────────────────────────────────────────────────────
            fv.ChargeState = result.ChargeState;
            fv.PrecursorMz = (float)result.PrecursorMz;
            fv.IsDecoy = result.IsDecoy;
            fv.PrecursorIndex = precursorIndex;

            return fv;
        }

        /// <summary>
        /// Looks up the observed RT at the apex time index from the XIC buffers.
        /// Uses the reference fragment (most data points) same as DiaTemporalScorer.
        /// </summary>
        private static float ComputeApexRt(
            ExtractionResult extractionResult,
            DiaLibraryQueryGenerator.PrecursorQueryGroup queryGroup,
            int apexTimeIndex)
        {
            int refFragIdx = -1;
            int maxPoints = 0;
            for (int f = 0; f < queryGroup.QueryCount; f++)
            {
                int qi = queryGroup.QueryOffset + f;
                int pts = extractionResult.Results[qi].DataPointCount;
                if (pts > maxPoints)
                {
                    maxPoints = pts;
                    refFragIdx = f;
                }
            }

            if (refFragIdx < 0 || maxPoints == 0 || apexTimeIndex < 0 || apexTimeIndex >= maxPoints)
                return float.NaN;

            int refQi = queryGroup.QueryOffset + refFragIdx;
            int rtOffset = extractionResult.Results[refQi].RtBufferOffset;
            return extractionResult.RtBuffer[rtOffset + apexTimeIndex];
        }

        /// <summary>Replaces NaN scores with 0 for classifier safety.</summary>
        private static float SafeScore(float score) =>
            float.IsNaN(score) ? 0f : score;
    }
}