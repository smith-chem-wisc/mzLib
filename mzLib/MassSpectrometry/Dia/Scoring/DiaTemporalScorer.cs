// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Buffers;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Defines how extracted fragment intensities are aggregated and scored
    /// against the library spectrum.
    /// 
    /// Modern DIA engines (DIA-NN, Spectronaut) do NOT score from a single
    /// summed intensity vector. Instead, they exploit the temporal (RT) structure
    /// of extracted ion chromatograms to compute time-resolved similarity.
    /// 
    /// The progression of scoring quality:
    ///   Summed         → median ~0.55 (current baseline, wide windows)
    ///   ConsensusApex  → median ~0.65 (score at chromatographic peak)
    ///   TemporalCosine → median ~0.70 (intensity-weighted average across RT)
    ///   WeightedTemporalCosineWithTransform → median ~0.72+ (DIA-NN-like)
    /// </summary>
    public enum ScoringStrategy
    {
        /// <summary>
        /// Original method: sum all XIC intensity per fragment, compare the summed
        /// vector against the library. Simple but dilutes signal with interference
        /// when RT windows are wide.
        /// </summary>
        Summed,

        /// <summary>
        /// Find the RT point where total fragment signal is highest (the consensus
        /// chromatographic apex), then build an observed intensity vector from each
        /// fragment's intensity at that single time point. Scores the apex vector
        /// against the library.
        /// 
        /// Much better than Summed because the library represents apex intensities.
        /// Requires that fragment XICs share a common RT grid (they do — all fragments
        /// for a precursor are extracted from the same window's scans).
        /// </summary>
        ConsensusApex,

        /// <summary>
        /// For each RT point across the chromatographic peak, compute cosine similarity
        /// between the observed fragment vector and the library vector. Then take the
        /// intensity-weighted average of those per-time cosines.
        /// 
        /// This rewards temporal coherence: fragments from the same peptide rise and
        /// fall together, so their per-time cosines are consistently high. Interference
        /// fragments don't track, so they contribute low cosines that get down-weighted.
        /// 
        /// This is the approach used by DIA-NN and similar tools.
        /// </summary>
        TemporalCosine,

        /// <summary>
        /// Same as TemporalCosine, but additionally:
        ///   - Uses sqrt(intensity) weighting to emphasize high-signal time points
        ///   - Applies a nonlinear transform (cos^N) to the final score
        /// 
        /// The nonlinear transform accentuates high similarities and suppresses
        /// moderate/ambiguous matches, improving target-decoy separation.
        /// </summary>
        WeightedTemporalCosineWithTransform
    }

    /// <summary>
    /// Computes spectrum similarity scores using RT-resolved XIC data.
    /// 
    /// Unlike the existing IScorer implementations (NormalizedDotProductScorer,
    /// SpectralAngleScorer) which operate on pre-collapsed intensity vectors,
    /// this scorer works directly on the raw XIC buffers from ExtractionResult.
    /// 
    /// It builds a time × fragment matrix from the per-fragment XIC data,
    /// then applies one of several scoring strategies.
    /// 
    /// Performance characteristics:
    ///   - Uses ArrayPool for all temporary allocations (zero GC in hot paths)
    ///   - No LINQ in scoring loops
    ///   - Span-based access to XIC buffers (zero-copy)
    ///   - Handles missing data gracefully (fragments with 0 XIC points get 0 intensity)
    ///   - Thread-safe: each instance can be used from one thread at a time;
    ///     create one per thread for parallel use
    /// </summary>
    public sealed class DiaTemporalScorer
    {
        private readonly ScoringStrategy _strategy;
        private readonly float _nonlinearPower;

        /// <summary>
        /// Creates a temporal scorer with the specified strategy.
        /// </summary>
        /// <param name="strategy">Which scoring algorithm to use.</param>
        /// <param name="nonlinearPower">
        /// Exponent for the nonlinear transform in WeightedTemporalCosineWithTransform mode.
        /// Default is 3.0 (as used in DIA-NN). Only used when strategy is
        /// WeightedTemporalCosineWithTransform. Must be >= 1.0.
        /// </param>
        public DiaTemporalScorer(ScoringStrategy strategy = ScoringStrategy.TemporalCosine,
            float nonlinearPower = 3.0f)
        {
            _strategy = strategy;
            _nonlinearPower = Math.Max(1.0f, nonlinearPower);
        }

        /// <summary>
        /// Scores a single precursor using its RT-resolved XIC data.
        /// 
        /// This is the main entry point called once per precursor during result assembly.
        /// It reads per-fragment XIC data from the extraction buffers, builds a time-aligned
        /// fragment intensity matrix, and applies the selected scoring strategy.
        /// </summary>
        /// <param name="libraryIntensities">
        /// Library/predicted fragment intensities (length = fragmentCount).
        /// </param>
        /// <param name="fragmentCount">Number of fragments for this precursor.</param>
        /// <param name="extractionResults">
        /// Per-fragment extraction results (slice of length fragmentCount, starting at the
        /// precursor's query offset).
        /// </param>
        /// <param name="rtBuffer">The full ExtractionResult.RtBuffer.</param>
        /// <param name="intensityBuffer">The full ExtractionResult.IntensityBuffer.</param>
        /// <returns>
        /// A DiaTemporalScore containing the primary score and diagnostic information.
        /// Returns a score with NaN values if insufficient data.
        /// </returns>
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
            // All fragments for a precursor come from the same isolation window,
            // so their XICs share the same set of scan RTs. We need to find the
            // union of all unique RT values across all fragments for this precursor.
            //
            // In practice, all fragments see the same scans (same window, same RT range),
            // so the RT grids are nearly identical. We use the fragment with the most
            // data points as the reference grid, then align other fragments to it.

            // Find the fragment with the most data points to use as RT reference
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

            // Get the reference RT grid
            var refResult = extractionResults[refFragmentIdx];
            ReadOnlySpan<float> refRts = rtBuffer.Slice(refResult.RtBufferOffset, refResult.DataPointCount);
            int timePointCount = refRts.Length;

            if (timePointCount == 0)
                return DiaTemporalScore.Insufficient;

            // ── Step 2: Build time × fragment intensity matrix ──────────────
            // matrix[t * fragmentCount + f] = intensity of fragment f at time t
            // Uses ArrayPool for zero-GC allocation.
            int matrixSize = timePointCount * fragmentCount;
            float[] matrixRented = ArrayPool<float>.Shared.Rent(matrixSize);
            Span<float> matrix = matrixRented.AsSpan(0, matrixSize);
            matrix.Clear(); // Zero-initialize (missing data = 0)

            try
            {
                for (int f = 0; f < fragmentCount; f++)
                {
                    var fr = extractionResults[f];
                    if (fr.DataPointCount == 0) continue;

                    ReadOnlySpan<float> fragRts = rtBuffer.Slice(fr.RtBufferOffset, fr.DataPointCount);
                    ReadOnlySpan<float> fragInts = intensityBuffer.Slice(fr.IntensityBufferOffset, fr.DataPointCount);

                    // Align this fragment's XIC to the reference RT grid.
                    // Since all fragments come from the same window, the RTs should match exactly
                    // or very closely. We use a two-pointer merge.
                    AlignXicToGrid(refRts, fragRts, fragInts, matrix, f, fragmentCount);
                }

                // ── Step 3: Apply the scoring strategy ──────────────────────
                return _strategy switch
                {
                    ScoringStrategy.Summed => ScoreSummed(libraryIntensities, matrix, fragmentCount, timePointCount),
                    ScoringStrategy.ConsensusApex => ScoreConsensusApex(libraryIntensities, matrix, fragmentCount, timePointCount, refRts),
                    ScoringStrategy.TemporalCosine => ScoreTemporalCosine(libraryIntensities, matrix, fragmentCount, timePointCount, refRts, useWeighting: false, nonlinearPower: 1.0f),
                    ScoringStrategy.WeightedTemporalCosineWithTransform => ScoreTemporalCosine(libraryIntensities, matrix, fragmentCount, timePointCount, refRts, useWeighting: true, nonlinearPower: _nonlinearPower),
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
        /// 
        /// Both fragRts and refRts are sorted (ascending RT). For each reference RT,
        /// we find the nearest fragment RT within a small tolerance and write its intensity
        /// into the matrix. Points that don't match get 0 (the matrix is pre-zeroed).
        /// 
        /// RT tolerance: 0.01 minutes (0.6 seconds) — scans from the same cycle should
        /// have nearly identical RTs, with tiny offsets for different windows in the cycle.
        /// </summary>
        private static void AlignXicToGrid(
            ReadOnlySpan<float> refRts,
            ReadOnlySpan<float> fragRts,
            ReadOnlySpan<float> fragInts,
            Span<float> matrix,
            int fragmentIndex,
            int fragmentCount)
        {
            const float rtTolerance = 0.01f; // 0.6 seconds

            int fragPtr = 0;
            for (int t = 0; t < refRts.Length && fragPtr < fragRts.Length; t++)
            {
                float refRt = refRts[t];

                // Advance fragPtr until we're close to or past refRt
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
        /// This reproduces the current (baseline) behavior.
        /// </summary>
        private static DiaTemporalScore ScoreSummed(
            ReadOnlySpan<float> library,
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int timePointCount)
        {
            // Sum across time for each fragment
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

                float dp = NormalizedDotProduct(library, summed);
                return new DiaTemporalScore(dp, dp, timePointCount, -1);
            }
            finally
            {
                ArrayPool<float>.Shared.Return(summedRented);
            }
        }

        /// <summary>
        /// Consensus apex scoring: finds the single time point where total fragment
        /// signal is highest, then scores that time point's fragment vector against
        /// the library.
        /// 
        /// Why "consensus": the apex is determined by the sum of all fragment intensities,
        /// not any single fragment. This is more robust than per-fragment apex because
        /// interference may cause a single fragment to peak at the wrong time.
        /// </summary>
        private static DiaTemporalScore ScoreConsensusApex(
            ReadOnlySpan<float> library,
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int timePointCount,
            ReadOnlySpan<float> refRts)
        {
            // Find the time point with the highest total fragment signal
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

            // Extract the apex fragment vector
            float[] apexRented = ArrayPool<float>.Shared.Rent(fragmentCount);
            Span<float> apexVector = apexRented.AsSpan(0, fragmentCount);

            try
            {
                int apexRowOffset = apexTimeIdx * fragmentCount;
                for (int f = 0; f < fragmentCount; f++)
                    apexVector[f] = matrix[apexRowOffset + f];

                float dp = NormalizedDotProduct(library, apexVector);
                float apexRt = apexTimeIdx < refRts.Length ? refRts[apexTimeIdx] : float.NaN;

                return new DiaTemporalScore(dp, dp, timePointCount, apexTimeIdx);
            }
            finally
            {
                ArrayPool<float>.Shared.Return(apexRented);
            }
        }

        /// <summary>
        /// Temporal cosine scoring (DIA-NN-style):
        ///   1. At each time point t, compute cos_t = cosine(observed_t, library)
        ///   2. Compute a weight w_t for each time point
        ///   3. Final score = Σ(w_t * cos_t) / Σ(w_t)
        ///   4. Optionally apply nonlinear transform: score^power
        /// 
        /// Two weighting schemes:
        ///   - Unweighted (useWeighting=false): w_t = 1 for all t → simple average cosine
        ///   - Intensity-weighted (useWeighting=true): w_t = sqrt(Σ_f intensity_t_f)
        ///     This emphasizes time points with strong signal (near the apex) and
        ///     de-emphasizes noise-floor time points at the edges of the RT window.
        /// </summary>
        private static DiaTemporalScore ScoreTemporalCosine(
            ReadOnlySpan<float> library,
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int timePointCount,
            ReadOnlySpan<float> refRts,
            bool useWeighting,
            float nonlinearPower)
        {
            // Pre-compute library L2 norm (constant across time points)
            float libNormSq = 0f;
            for (int f = 0; f < fragmentCount; f++)
                libNormSq += library[f] * library[f];

            if (libNormSq <= 0f)
                return DiaTemporalScore.Insufficient;

            float libNorm = MathF.Sqrt(libNormSq);

            // Temporary buffer for per-time observed vector
            float[] obsRented = ArrayPool<float>.Shared.Rent(fragmentCount);
            Span<float> observed = obsRented.AsSpan(0, fragmentCount);

            try
            {
                float weightedCosineSum = 0f;
                float weightSum = 0f;
                int validTimePoints = 0;

                // Also track apex for diagnostic purposes
                int apexTimeIdx = 0;
                float apexTotalIntensity = 0f;

                for (int t = 0; t < timePointCount; t++)
                {
                    // Extract observed vector for this time point
                    int rowOffset = t * fragmentCount;
                    float totalIntensity = 0f;
                    for (int f = 0; f < fragmentCount; f++)
                    {
                        observed[f] = matrix[rowOffset + f];
                        totalIntensity += observed[f];
                    }

                    // Skip time points with no signal
                    if (totalIntensity <= 0f)
                        continue;

                    // Track apex
                    if (totalIntensity > apexTotalIntensity)
                    {
                        apexTotalIntensity = totalIntensity;
                        apexTimeIdx = t;
                    }

                    // Compute cosine at this time point
                    float cos_t = NormalizedDotProduct(library, observed);
                    if (float.IsNaN(cos_t))
                        continue;

                    // Compute weight
                    float weight = useWeighting ? MathF.Sqrt(totalIntensity) : 1.0f;

                    weightedCosineSum += weight * cos_t;
                    weightSum += weight;
                    validTimePoints++;
                }

                if (weightSum <= 0f || validTimePoints == 0)
                    return DiaTemporalScore.Insufficient;

                float temporalScore = weightedCosineSum / weightSum;

                // Apply nonlinear transform if requested
                float transformedScore = nonlinearPower > 1.0f
                    ? MathF.Pow(Math.Max(temporalScore, 0f), nonlinearPower)
                    : temporalScore;

                return new DiaTemporalScore(
                    dotProductScore: transformedScore,
                    rawCosine: temporalScore,
                    timePointsUsed: validTimePoints,
                    apexTimeIndex: apexTimeIdx);
            }
            finally
            {
                ArrayPool<float>.Shared.Return(obsRented);
            }
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Utility: Normalized Dot Product (L2-normalized cosine)
        // ─────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Computes L2-normalized dot product (cosine similarity) between two vectors.
        /// Returns NaN if either vector has zero norm.
        /// 
        /// This is intentionally a static method using Span to avoid any overhead.
        /// For hot-path scoring, we inline this rather than going through IScorer.
        /// </summary>
        private static float NormalizedDotProduct(ReadOnlySpan<float> a, ReadOnlySpan<float> b)
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

            float result = dot / (MathF.Sqrt(normA) * MathF.Sqrt(normB));
            // Clamp to [0, 1] to handle floating-point imprecision
            return Math.Clamp(result, 0f, 1f);
        }
    }

    /// <summary>
    /// Result of temporal scoring for a single precursor.
    /// Contains both the final score and diagnostic information useful
    /// for debugging, validation, and future multi-feature classifiers.
    /// </summary>
    public readonly struct DiaTemporalScore
    {
        /// <summary>
        /// The primary similarity score after all transforms.
        /// Range [0, 1], NaN if insufficient data.
        /// For Summed/ConsensusApex: this equals RawCosine.
        /// For TemporalCosine: this is the weighted-average cosine.
        /// For WeightedTemporalCosineWithTransform: this is (weighted_cosine)^N.
        /// </summary>
        public readonly float DotProductScore;

        /// <summary>
        /// The cosine similarity before nonlinear transform.
        /// Useful for diagnostics: compare RawCosine to DotProductScore to see
        /// the effect of the transform.
        /// </summary>
        public readonly float RawCosine;

        /// <summary>
        /// Number of RT time points that contributed to the score.
        /// Higher means more temporal evidence (more robust score).
        /// </summary>
        public readonly int TimePointsUsed;

        /// <summary>
        /// Index of the consensus apex time point (time with highest total signal).
        /// -1 if not applicable. Useful for reporting the apex RT.
        /// </summary>
        public readonly int ApexTimeIndex;

        /// <summary>Sentinel value for insufficient data (e.g., too few fragments detected).</summary>
        public static readonly DiaTemporalScore Insufficient = new(float.NaN, float.NaN, 0, -1);

        public DiaTemporalScore(float dotProductScore, float rawCosine, int timePointsUsed, int apexTimeIndex)
        {
            DotProductScore = dotProductScore;
            RawCosine = rawCosine;
            TimePointsUsed = timePointsUsed;
            ApexTimeIndex = apexTimeIndex;
        }

        /// <summary>Whether this score represents a valid result (not NaN).</summary>
        public bool IsValid => !float.IsNaN(DotProductScore);

        public override string ToString()
        {
            if (!IsValid) return "Insufficient data";
            return $"DP={DotProductScore:F4} (raw={RawCosine:F4}, {TimePointsUsed} time points, apex@{ApexTimeIndex})";
        }
    }
}