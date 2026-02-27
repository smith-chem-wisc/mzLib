// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Buffers;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Defines how extracted fragment intensities are aggregated and scored
    /// against the library spectrum.
    /// </summary>
    public enum ScoringStrategy
    {
        /// <summary>
        /// Original method: sum all XIC intensity per fragment, compare the summed
        /// vector against the library.
        /// </summary>
        Summed,

        /// <summary>
        /// Find the RT point where total fragment signal is highest (the consensus
        /// chromatographic apex), then build an observed intensity vector from each
        /// fragment's intensity at that single time point.
        /// </summary>
        ConsensusApex,

        /// <summary>
        /// For each RT point across the chromatographic peak, compute cosine similarity
        /// between the observed fragment vector and the library vector. Then take the
        /// intensity-weighted average of those per-time cosines.
        /// 
        /// At each time point, only fragments with nonzero signal participate in
        /// the cosine calculation. This prevents missing fragments (due to no m/z
        /// match in a given scan) from pulling down the score.
        /// </summary>
        TemporalCosine,

        /// <summary>
        /// Same as TemporalCosine, but additionally:
        ///   - Uses sqrt(intensity) weighting to emphasize high-signal time points
        ///   - Applies a nonlinear transform (cos^N) to the final score
        /// </summary>
        WeightedTemporalCosineWithTransform
    }

    /// <summary>
    /// Computes spectrum similarity scores using RT-resolved XIC data.
    /// 
    /// KEY DESIGN DECISION: In DIA extraction, each fragment independently searches
    /// for its target m/z in each scan. A fragment only appears in the XIC at a given
    /// RT if a matching peak was found. Therefore, different fragments may have data
    /// at different subsets of time points. This scorer handles this by:
    ///   - Building a time × fragment matrix aligned to a common RT grid
    ///   - At each time point, computing cosine using only fragments WITH signal
    ///   - Requiring a minimum number of active fragments per time point
    ///   - Weighting time points by signal strength
    /// 
    /// Thread-safe: create one instance per thread for parallel use.
    /// </summary>
    public sealed class DiaTemporalScorer
    {
        private readonly ScoringStrategy _strategy;
        private readonly float _nonlinearPower;
        private readonly int _minActiveFragments;

        /// <summary>
        /// Creates a temporal scorer with the specified strategy.
        /// </summary>
        /// <param name="strategy">Which scoring algorithm to use.</param>
        /// <param name="nonlinearPower">
        /// Exponent for WeightedTemporalCosineWithTransform. Default 3.0.
        /// </param>
        /// <param name="minActiveFragments">
        /// Minimum number of fragments that must have nonzero intensity at a time point
        /// for that point to participate in temporal scoring. Default 3.
        /// Time points with fewer active fragments are too noisy to produce meaningful
        /// cosine values and are skipped.
        /// </param>
        public DiaTemporalScorer(ScoringStrategy strategy = ScoringStrategy.TemporalCosine,
            float nonlinearPower = 3.0f,
            int minActiveFragments = 3)
        {
            _strategy = strategy;
            _nonlinearPower = Math.Max(1.0f, nonlinearPower);
            _minActiveFragments = Math.Max(2, minActiveFragments);
        }

        /// <summary>
        /// Scores a single precursor using its RT-resolved XIC data.
        /// </summary>
        public DiaTemporalScore ScorePrecursor(
            ReadOnlySpan<float> libraryIntensities,
            int fragmentCount,
            ReadOnlySpan<FragmentResult> extractionResults,
            ReadOnlySpan<float> rtBuffer,
            ReadOnlySpan<float> intensityBuffer)
        {
            if (fragmentCount < 2)
                return DiaTemporalScore.Insufficient;

            // ── Step 1: Discover the common RT grid ─────────────────────────
            int refFragmentIdx = -1;
            int maxPoints = 0;
            for (int f = 0; f < fragmentCount; f++)
            {
                if (extractionResults[f].DataPointCount > maxPoints)
                {
                    maxPoints = extractionResults[f].DataPointCount;
                    refFragmentIdx = f;
                }
            }

            if (maxPoints == 0)
                return DiaTemporalScore.Insufficient;

            var refResult = extractionResults[refFragmentIdx];
            ReadOnlySpan<float> refRts = rtBuffer.Slice(refResult.RtBufferOffset, refResult.DataPointCount);
            int timePointCount = refRts.Length;

            if (timePointCount == 0)
                return DiaTemporalScore.Insufficient;

            // ── Step 2: Build time × fragment intensity matrix ──────────────
            int matrixSize = timePointCount * fragmentCount;
            float[] matrixRented = ArrayPool<float>.Shared.Rent(matrixSize);
            Span<float> matrix = matrixRented.AsSpan(0, matrixSize);
            matrix.Clear();

            try
            {
                for (int f = 0; f < fragmentCount; f++)
                {
                    var fr = extractionResults[f];
                    if (fr.DataPointCount == 0) continue;

                    ReadOnlySpan<float> fragRts = rtBuffer.Slice(fr.RtBufferOffset, fr.DataPointCount);
                    ReadOnlySpan<float> fragInts = intensityBuffer.Slice(fr.IntensityBufferOffset, fr.DataPointCount);

                    AlignXicToGrid(refRts, fragRts, fragInts, matrix, f, fragmentCount);
                }

                return _strategy switch
                {
                    ScoringStrategy.Summed => ScoreSummed(libraryIntensities, matrix, fragmentCount, timePointCount),
                    ScoringStrategy.ConsensusApex => ScoreConsensusApex(libraryIntensities, matrix, fragmentCount, timePointCount),
                    ScoringStrategy.TemporalCosine => ScoreTemporalCosine(libraryIntensities, matrix, fragmentCount, timePointCount, useWeighting: false, nonlinearPower: 1.0f),
                    ScoringStrategy.WeightedTemporalCosineWithTransform => ScoreTemporalCosine(libraryIntensities, matrix, fragmentCount, timePointCount, useWeighting: true, nonlinearPower: _nonlinearPower),
                    _ => DiaTemporalScore.Insufficient
                };
            }
            finally
            {
                ArrayPool<float>.Shared.Return(matrixRented);
            }
        }

        /// <summary>
        /// Aligns a fragment's XIC data points to the reference RT grid using two-pointer merge.
        /// RT tolerance: 0.01 minutes (0.6 seconds).
        /// </summary>
        private static void AlignXicToGrid(
            ReadOnlySpan<float> refRts,
            ReadOnlySpan<float> fragRts,
            ReadOnlySpan<float> fragInts,
            Span<float> matrix,
            int fragmentIndex,
            int fragmentCount)
        {
            const float rtTolerance = 0.01f;

            int fragPtr = 0;
            for (int t = 0; t < refRts.Length && fragPtr < fragRts.Length; t++)
            {
                float refRt = refRts[t];

                while (fragPtr < fragRts.Length && fragRts[fragPtr] < refRt - rtTolerance)
                    fragPtr++;

                if (fragPtr < fragRts.Length && MathF.Abs(fragRts[fragPtr] - refRt) <= rtTolerance)
                {
                    matrix[t * fragmentCount + fragmentIndex] = fragInts[fragPtr];
                    fragPtr++;
                }
            }
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Scoring Strategy Implementations
        // ─────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Summed scoring: collapses all time points by summing, then computes cosine.
        /// </summary>
        private static DiaTemporalScore ScoreSummed(
            ReadOnlySpan<float> library,
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int timePointCount)
        {
            float[] summedRented = ArrayPool<float>.Shared.Rent(fragmentCount);
            Span<float> summed = summedRented.AsSpan(0, fragmentCount);
            summed.Clear();

            try
            {
                for (int t = 0; t < timePointCount; t++)
                {
                    int rowOffset = t * fragmentCount;
                    for (int f = 0; f < fragmentCount; f++)
                        summed[f] += matrix[rowOffset + f];
                }

                float dp = CosineAllFragments(library, summed);
                return new DiaTemporalScore(dp, dp, timePointCount, -1);
            }
            finally
            {
                ArrayPool<float>.Shared.Return(summedRented);
            }
        }

        /// <summary>
        /// Consensus apex scoring: finds the single time point where total fragment
        /// signal is highest, then scores that time point's fragment vector.
        /// Only includes fragments with nonzero intensity at the apex in the cosine.
        /// </summary>
        private static DiaTemporalScore ScoreConsensusApex(
            ReadOnlySpan<float> library,
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int timePointCount)
        {
            int apexTimeIdx = 0;
            float apexTotalIntensity = 0f;

            for (int t = 0; t < timePointCount; t++)
            {
                float total = 0f;
                int rowOffset = t * fragmentCount;
                for (int f = 0; f < fragmentCount; f++)
                    total += matrix[rowOffset + f];

                if (total > apexTotalIntensity)
                {
                    apexTotalIntensity = total;
                    apexTimeIdx = t;
                }
            }

            if (apexTotalIntensity <= 0f)
                return DiaTemporalScore.Insufficient;

            int apexRowOffset = apexTimeIdx * fragmentCount;
            float dp = CosineActiveFragments(library, matrix, apexRowOffset, fragmentCount);

            return new DiaTemporalScore(dp, dp, timePointCount, apexTimeIdx);
        }

        /// <summary>
        /// Temporal cosine scoring:
        ///   1. At each time point, compute cosine using ONLY fragments with nonzero intensity
        ///   2. Require at least _minActiveFragments per time point
        ///   3. Weight by total signal (or uniform)
        ///   4. Optionally apply nonlinear transform
        /// 
        /// The critical fix: zero-intensity entries (fragments not found in a given scan)
        /// are EXCLUDED from the cosine calculation, not treated as 0.
        /// </summary>
        private DiaTemporalScore ScoreTemporalCosine(
            ReadOnlySpan<float> library,
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int timePointCount,
            bool useWeighting,
            float nonlinearPower)
        {
            float weightedCosineSum = 0f;
            float weightSum = 0f;
            int validTimePoints = 0;
            int apexTimeIdx = 0;
            float apexTotalIntensity = 0f;

            for (int t = 0; t < timePointCount; t++)
            {
                int rowOffset = t * fragmentCount;

                // Count active fragments and total intensity
                int activeCount = 0;
                float totalIntensity = 0f;
                for (int f = 0; f < fragmentCount; f++)
                {
                    if (matrix[rowOffset + f] > 0f)
                    {
                        activeCount++;
                        totalIntensity += matrix[rowOffset + f];
                    }
                }

                // Skip time points with too few active fragments
                if (activeCount < _minActiveFragments || totalIntensity <= 0f)
                    continue;

                // Track apex
                if (totalIntensity > apexTotalIntensity)
                {
                    apexTotalIntensity = totalIntensity;
                    apexTimeIdx = t;
                }

                // Compute cosine using only active fragments
                float cos_t = CosineActiveFragments(library, matrix, rowOffset, fragmentCount);
                if (float.IsNaN(cos_t))
                    continue;

                float weight = useWeighting ? MathF.Sqrt(totalIntensity) : 1.0f;

                weightedCosineSum += weight * cos_t;
                weightSum += weight;
                validTimePoints++;
            }

            if (weightSum <= 0f || validTimePoints == 0)
                return DiaTemporalScore.Insufficient;

            float temporalScore = weightedCosineSum / weightSum;

            float transformedScore = nonlinearPower > 1.0f
                ? MathF.Pow(Math.Max(temporalScore, 0f), nonlinearPower)
                : temporalScore;

            return new DiaTemporalScore(
                dotProductScore: transformedScore,
                rawCosine: temporalScore,
                timePointsUsed: validTimePoints,
                apexTimeIndex: apexTimeIdx);
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Cosine Similarity Helpers
        // ─────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Standard cosine between two full vectors (for summed scoring).
        /// </summary>
        private static float CosineAllFragments(ReadOnlySpan<float> a, ReadOnlySpan<float> b)
        {
            int len = Math.Min(a.Length, b.Length);
            float dot = 0f, normA = 0f, normB = 0f;

            for (int i = 0; i < len; i++)
            {
                dot += a[i] * b[i];
                normA += a[i] * a[i];
                normB += b[i] * b[i];
            }

            if (normA <= 0f || normB <= 0f)
                return float.NaN;

            return Math.Clamp(dot / (MathF.Sqrt(normA) * MathF.Sqrt(normB)), 0f, 1f);
        }

        /// <summary>
        /// Cosine between library and observed, using ONLY fragments where observed > 0.
        /// 
        /// This is the key to correct temporal scoring: at a given time point, some
        /// fragments may not have been found (observed = 0) because no peak matched
        /// their m/z in that scan. Including these zeros in the cosine would falsely
        /// reduce the score. Instead, we compute cosine only over the subset of
        /// fragments that actually have signal.
        /// </summary>
        private static float CosineActiveFragments(
            ReadOnlySpan<float> library,
            ReadOnlySpan<float> matrix,
            int rowOffset,
            int fragmentCount)
        {
            float dot = 0f, normLib = 0f, normObs = 0f;

            for (int f = 0; f < fragmentCount; f++)
            {
                float obs = matrix[rowOffset + f];
                if (obs <= 0f) continue;

                float lib = f < library.Length ? library[f] : 0f;
                dot += lib * obs;
                normLib += lib * lib;
                normObs += obs * obs;
            }

            if (normLib <= 0f || normObs <= 0f)
                return float.NaN;

            return Math.Clamp(dot / (MathF.Sqrt(normLib) * MathF.Sqrt(normObs)), 0f, 1f);
        }
    }

    /// <summary>
    /// Result of temporal scoring for a single precursor.
    /// </summary>
    public readonly struct DiaTemporalScore
    {
        public readonly float DotProductScore;
        public readonly float RawCosine;
        public readonly int TimePointsUsed;
        public readonly int ApexTimeIndex;

        public static readonly DiaTemporalScore Insufficient = new(float.NaN, float.NaN, 0, -1);

        public DiaTemporalScore(float dotProductScore, float rawCosine, int timePointsUsed, int apexTimeIndex)
        {
            DotProductScore = dotProductScore;
            RawCosine = rawCosine;
            TimePointsUsed = timePointsUsed;
            ApexTimeIndex = apexTimeIndex;
        }

        public bool IsValid => !float.IsNaN(DotProductScore);

        public override string ToString()
        {
            if (!IsValid) return "Insufficient data";
            return $"DP={DotProductScore:F4} (raw={RawCosine:F4}, {TimePointsUsed} time points, apex@{ApexTimeIndex})";
        }
    }
}