// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
// Location: MassSpectrometry/Dia/Calibration/RobustApexDetector.cs

using System;
using System.Buffers;
using System.Runtime.CompilerServices;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Robust apex detection for DIA precursors that resists interference contamination.
    /// 
    /// Problem: The standard apex detection (max total fragment intensity across time points)
    /// is highly susceptible to interference from co-eluting peptides. In a wide RT window
    /// (±5.0 min), the "apex" frequently corresponds to an interfering peptide's peak rather
    /// than the target's true elution time. This corrupts RT calibration anchors.
    /// 
    /// Solution: Use the "best fragment" — the fragment whose XIC correlates most strongly
    /// with all other fragments (identified by DiaFeatureCalculator.ComputeBestFragmentFeatures).
    /// This fragment is least affected by interference. By finding the apex of just this one
    /// fragment's smoothed XIC, we get a much more reliable estimate of the target's true
    /// elution time.
    /// 
    /// For ambiguous cases (multiple peaks in the best fragment's XIC), the multi-candidate
    /// method identifies all local maxima and scores each by a composite of intensity and
    /// library spectral agreement, selecting the candidate most likely to be the true target.
    /// 
    /// Performance: O(T) per precursor for single-apex, O(T × F) for multi-candidate scoring
    /// where T = time points, F = fragments. No heap allocations (uses stackalloc/ArrayPool).
    /// Thread-safe: stateless, operates only on input spans.
    /// </summary>
    public static class RobustApexDetector
    {
        /// <summary>
        /// Gaussian smoothing sigma in scan units. σ=1 means the kernel extends ±3 scans
        /// with significant weight, providing good noise suppression without over-smoothing
        /// typical chromatographic peaks (which span 5-15 scans).
        /// </summary>
        private const float SmoothingSigma = 1.0f;

        /// <summary>
        /// Half-width of the Gaussian kernel in scan units.
        /// 3σ captures >99% of the kernel's mass.
        /// </summary>
        private const int KernelHalfWidth = 3;

        /// <summary>
        /// Minimum fraction of the global max in the smoothed XIC for a local maximum
        /// to be considered a real apex candidate (not a noise spike).
        /// </summary>
        private const float MinPeakFraction = 0.10f;

        // ─────────────────────────────────────────────────────────────────────
        //  Task 2: Robust Single-Apex Detection
        // ─────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Identifies the apex scan using only the best (least-interfered) fragment's XIC,
        /// rather than the sum across all fragments. The best fragment is the one whose
        /// XIC correlates most strongly with all other fragments (from ComputeBestFragmentFeatures).
        /// 
        /// Falls back to TIC-based apex if bestFragmentIndex is unavailable (&lt;0) or
        /// if fewer than 3 time points have signal in the best fragment.
        /// 
        /// Algorithm:
        ///   1. Extract the best fragment's column from the time×fragment matrix
        ///   2. Apply Gaussian smoothing (σ = 1 scan) to reduce noise spikes
        ///   3. Find the time point with maximum smoothed intensity
        ///   4. If best fragment is unavailable, fall back to max-TIC apex
        /// </summary>
        /// <param name="matrix">
        /// Row-major time × fragment intensity matrix.
        /// Layout: matrix[t * fragmentCount + f] = intensity of fragment f at time point t.
        /// </param>
        /// <param name="fragmentCount">Number of fragments (columns).</param>
        /// <param name="timePointCount">Number of time points (rows).</param>
        /// <param name="bestFragmentIndex">
        /// Index of the best fragment from DiaFeatureCalculator/DiaBestFragmentHelper.
        /// Pass -1 if not available (will use TIC fallback).
        /// </param>
        /// <returns>Zero-based time-point index of the robust apex.</returns>
        public static int FindRobustApex(
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int timePointCount,
            int bestFragmentIndex)
        {
            if (timePointCount <= 0) return 0;
            if (timePointCount == 1) return 0;

            // If best fragment is available and has signal, use it
            if (bestFragmentIndex >= 0 && bestFragmentIndex < fragmentCount)
            {
                int nonzero = CountNonzeroInColumn(matrix, bestFragmentIndex, fragmentCount, timePointCount);
                if (nonzero >= 3)
                {
                    return FindSmoothedApexForFragment(matrix, bestFragmentIndex, fragmentCount, timePointCount);
                }
            }

            // Fallback: TIC-based apex (sum all fragments per time point, find max)
            return FindTicApex(matrix, fragmentCount, timePointCount);
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Task 3: Multi-Apex Candidate Scoring
        // ─────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Returns up to maxCandidates apex candidates from the best fragment's smoothed XIC,
        /// each scored by a composite of (smoothed intensity at candidate) × (cosine similarity
        /// to library intensities at that candidate scan).
        /// 
        /// The caller should pick the candidate with the highest CompositeScore as the final apex.
        /// This handles the common case where target and interferer produce two distinct peaks:
        /// the target peak will have better library agreement even if the interferer is brighter.
        /// 
        /// Algorithm:
        ///   1. Extract and Gaussian-smooth the best fragment's XIC
        ///   2. Find all local maxima above MinPeakFraction of global max
        ///   3. For each local maximum, compute cosine similarity between the observed
        ///      fragment intensities at that time point and the library intensities
        ///   4. Score each candidate as: smoothedIntensity × cosineSimilarity
        ///   5. Return top-N candidates sorted by composite score (descending)
        /// 
        /// If bestFragmentIndex is unavailable, uses TIC profile instead.
        /// If fewer candidates than maxCandidates are found, returns only those found.
        /// Returns an empty array if no valid candidates exist.
        /// </summary>
        /// <param name="matrix">Row-major [timePoints × fragments] intensity matrix.</param>
        /// <param name="libraryIntensities">
        /// Library/predicted intensities per fragment (length = fragmentCount).
        /// Used for cosine scoring at each candidate time point.
        /// </param>
        /// <param name="fragmentCount">Number of fragments (columns).</param>
        /// <param name="timePointCount">Number of time points (rows).</param>
        /// <param name="bestFragmentIndex">
        /// Index of the best fragment. Pass -1 to use TIC profile instead.
        /// </param>
        /// <param name="maxCandidates">Maximum number of candidates to return (default 3).</param>
        /// <returns>
        /// Array of (ScanIndex, CompositeScore) tuples sorted by CompositeScore descending.
        /// ScanIndex is the zero-based time-point index into the matrix.
        /// </returns>
        public static (int ScanIndex, float CompositeScore)[] FindApexCandidates(
            ReadOnlySpan<float> matrix,
            ReadOnlySpan<float> libraryIntensities,
            int fragmentCount,
            int timePointCount,
            int bestFragmentIndex,
            int maxCandidates = 3)
        {
            if (timePointCount < 3 || fragmentCount < 1)
                return Array.Empty<(int, float)>();

            // Step 1: Build the smoothed profile (best fragment or TIC)
            float[] smoothedRented = ArrayPool<float>.Shared.Rent(timePointCount);
            Span<float> smoothed = smoothedRented.AsSpan(0, timePointCount);

            try
            {
                bool useBestFrag = bestFragmentIndex >= 0
                    && bestFragmentIndex < fragmentCount
                    && CountNonzeroInColumn(matrix, bestFragmentIndex, fragmentCount, timePointCount) >= 3;

                if (useBestFrag)
                {
                    GaussianSmoothColumn(matrix, bestFragmentIndex, fragmentCount, timePointCount, smoothed);
                }
                else
                {
                    // Build TIC profile then smooth it
                    float[] ticRented = ArrayPool<float>.Shared.Rent(timePointCount);
                    Span<float> tic = ticRented.AsSpan(0, timePointCount);
                    try
                    {
                        BuildTicProfile(matrix, fragmentCount, timePointCount, tic);
                        GaussianSmooth(tic, timePointCount, smoothed);
                    }
                    finally
                    {
                        ArrayPool<float>.Shared.Return(ticRented);
                    }
                }

                // Step 2: Find global max for thresholding
                float globalMax = 0f;
                for (int t = 0; t < timePointCount; t++)
                    if (smoothed[t] > globalMax) globalMax = smoothed[t];

                if (globalMax <= 0f)
                    return Array.Empty<(int, float)>();

                float minPeakHeight = globalMax * MinPeakFraction;

                // Step 3: Find all local maxima above threshold
                // Use stackalloc for candidate storage (limited to a reasonable max)
                const int maxLocalMaxima = 16;
                Span<int> localMaxIndices = stackalloc int[maxLocalMaxima];
                Span<float> localMaxValues = stackalloc float[maxLocalMaxima];
                int nLocalMax = 0;

                for (int t = 1; t < timePointCount - 1 && nLocalMax < maxLocalMaxima; t++)
                {
                    if (smoothed[t] >= minPeakHeight
                        && smoothed[t] >= smoothed[t - 1]
                        && smoothed[t] >= smoothed[t + 1]
                        && (smoothed[t] > smoothed[t - 1] || smoothed[t] > smoothed[t + 1]))
                    {
                        localMaxIndices[nLocalMax] = t;
                        localMaxValues[nLocalMax] = smoothed[t];
                        nLocalMax++;
                    }
                }

                // Also check endpoints if they're above threshold
                if (nLocalMax < maxLocalMaxima && timePointCount >= 2)
                {
                    if (smoothed[0] >= minPeakHeight && smoothed[0] >= smoothed[1])
                    {
                        localMaxIndices[nLocalMax] = 0;
                        localMaxValues[nLocalMax] = smoothed[0];
                        nLocalMax++;
                    }
                }
                if (nLocalMax < maxLocalMaxima && timePointCount >= 2)
                {
                    int last = timePointCount - 1;
                    if (smoothed[last] >= minPeakHeight && smoothed[last] >= smoothed[last - 1])
                    {
                        localMaxIndices[nLocalMax] = last;
                        localMaxValues[nLocalMax] = smoothed[last];
                        nLocalMax++;
                    }
                }

                if (nLocalMax == 0)
                {
                    // No local max found — return the global max as sole candidate
                    int bestT = 0;
                    for (int t = 1; t < timePointCount; t++)
                        if (smoothed[t] > smoothed[bestT]) bestT = t;

                    float cosine = ComputeCosineAtTimePoint(matrix, libraryIntensities, fragmentCount, bestT);
                    float score = smoothed[bestT] * Math.Max(cosine, 0f);
                    return new[] { (bestT, score) };
                }

                // Step 4: Score each candidate with composite = smoothedIntensity × cosine
                Span<float> compositeScores = stackalloc float[nLocalMax];
                for (int c = 0; c < nLocalMax; c++)
                {
                    int t = localMaxIndices[c];
                    float cosine = ComputeCosineAtTimePoint(matrix, libraryIntensities, fragmentCount, t);
                    // Clamp cosine to [0,1] — negative cosine means terrible match
                    compositeScores[c] = localMaxValues[c] * Math.Max(cosine, 0f);
                }

                // Step 5: Sort candidates by composite score (descending), take top-N
                // Simple selection-sort for small N
                int resultCount = Math.Min(nLocalMax, maxCandidates);
                var result = new (int ScanIndex, float CompositeScore)[resultCount];

                // Create index array for sorting
                Span<int> order = stackalloc int[nLocalMax];
                for (int i = 0; i < nLocalMax; i++) order[i] = i;

                // Partial selection sort: find top resultCount
                for (int i = 0; i < resultCount; i++)
                {
                    int bestIdx = i;
                    float bestScore = compositeScores[order[i]];
                    for (int j = i + 1; j < nLocalMax; j++)
                    {
                        if (compositeScores[order[j]] > bestScore)
                        {
                            bestScore = compositeScores[order[j]];
                            bestIdx = j;
                        }
                    }
                    // Swap
                    if (bestIdx != i)
                    {
                        int tmp = order[i];
                        order[i] = order[bestIdx];
                        order[bestIdx] = tmp;
                    }

                    int candIdx = order[i];
                    result[i] = (localMaxIndices[candIdx], compositeScores[candIdx]);
                }

                return result;
            }
            finally
            {
                ArrayPool<float>.Shared.Return(smoothedRented);
            }
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Internal Helpers
        // ─────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Finds the apex of a single fragment's Gaussian-smoothed XIC.
        /// Returns the time-point index with the highest smoothed intensity.
        /// </summary>
        private static int FindSmoothedApexForFragment(
            ReadOnlySpan<float> matrix,
            int fragmentIndex,
            int fragmentCount,
            int timePointCount)
        {
            // Rent buffer for smoothed values
            float[] smoothedRented = ArrayPool<float>.Shared.Rent(timePointCount);
            Span<float> smoothed = smoothedRented.AsSpan(0, timePointCount);

            try
            {
                GaussianSmoothColumn(matrix, fragmentIndex, fragmentCount, timePointCount, smoothed);

                int bestT = 0;
                float bestVal = smoothed[0];
                for (int t = 1; t < timePointCount; t++)
                {
                    if (smoothed[t] > bestVal)
                    {
                        bestVal = smoothed[t];
                        bestT = t;
                    }
                }
                return bestT;
            }
            finally
            {
                ArrayPool<float>.Shared.Return(smoothedRented);
            }
        }

        /// <summary>
        /// TIC-based apex: sum all fragment intensities per time point, return the max.
        /// This is the fallback when best-fragment is unavailable.
        /// </summary>
        private static int FindTicApex(
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int timePointCount)
        {
            int bestT = 0;
            float bestTic = 0f;

            for (int t = 0; t < timePointCount; t++)
            {
                float tic = 0f;
                int rowOffset = t * fragmentCount;
                for (int f = 0; f < fragmentCount; f++)
                    tic += matrix[rowOffset + f];

                if (tic > bestTic)
                {
                    bestTic = tic;
                    bestT = t;
                }
            }
            return bestT;
        }

        /// <summary>
        /// Builds the TIC (total ion current) profile by summing all fragments at each time point.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void BuildTicProfile(
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int timePointCount,
            Span<float> tic)
        {
            for (int t = 0; t < timePointCount; t++)
            {
                float sum = 0f;
                int rowOffset = t * fragmentCount;
                for (int f = 0; f < fragmentCount; f++)
                    sum += matrix[rowOffset + f];
                tic[t] = sum;
            }
        }

        /// <summary>
        /// Gaussian smoothing of a single column (fragment) from the row-major matrix.
        /// The kernel has σ = SmoothingSigma and extends ±KernelHalfWidth scans.
        /// Output is written to the provided span.
        /// </summary>
        private static void GaussianSmoothColumn(
            ReadOnlySpan<float> matrix,
            int fragmentIndex,
            int fragmentCount,
            int timePointCount,
            Span<float> output)
        {
            // Precompute Gaussian kernel weights
            Span<float> kernel = stackalloc float[2 * KernelHalfWidth + 1];
            float kernelSum = 0f;
            float twoSigmaSq = 2.0f * SmoothingSigma * SmoothingSigma;

            for (int k = -KernelHalfWidth; k <= KernelHalfWidth; k++)
            {
                float w = MathF.Exp(-(k * k) / twoSigmaSq);
                kernel[k + KernelHalfWidth] = w;
                kernelSum += w;
            }

            // Normalize kernel
            for (int k = 0; k < kernel.Length; k++)
                kernel[k] /= kernelSum;

            // Apply convolution
            for (int t = 0; t < timePointCount; t++)
            {
                float val = 0f;
                float wsum = 0f;

                for (int k = -KernelHalfWidth; k <= KernelHalfWidth; k++)
                {
                    int idx = t + k;
                    if (idx < 0 || idx >= timePointCount) continue;

                    float w = kernel[k + KernelHalfWidth];
                    val += matrix[idx * fragmentCount + fragmentIndex] * w;
                    wsum += w;
                }

                output[t] = wsum > 0f ? val / wsum : 0f;
            }
        }

        /// <summary>
        /// Gaussian smoothing of a 1D signal (e.g., pre-built TIC profile).
        /// Uses the same kernel parameters as the column smoother.
        /// </summary>
        private static void GaussianSmooth(
            ReadOnlySpan<float> input,
            int length,
            Span<float> output)
        {
            Span<float> kernel = stackalloc float[2 * KernelHalfWidth + 1];
            float kernelSum = 0f;
            float twoSigmaSq = 2.0f * SmoothingSigma * SmoothingSigma;

            for (int k = -KernelHalfWidth; k <= KernelHalfWidth; k++)
            {
                float w = MathF.Exp(-(k * k) / twoSigmaSq);
                kernel[k + KernelHalfWidth] = w;
                kernelSum += w;
            }

            for (int k = 0; k < kernel.Length; k++)
                kernel[k] /= kernelSum;

            for (int t = 0; t < length; t++)
            {
                float val = 0f;
                float wsum = 0f;

                for (int k = -KernelHalfWidth; k <= KernelHalfWidth; k++)
                {
                    int idx = t + k;
                    if (idx < 0 || idx >= length) continue;

                    float w = kernel[k + KernelHalfWidth];
                    val += input[idx] * w;
                    wsum += w;
                }

                output[t] = wsum > 0f ? val / wsum : 0f;
            }
        }

        /// <summary>
        /// Computes cosine similarity between library intensities and the observed
        /// fragment intensities at a specific time point in the matrix.
        /// 
        /// cosine = dot(lib, obs) / (||lib|| × ||obs||)
        /// 
        /// Only includes fragments where both library and observed intensities are positive.
        /// Returns NaN if fewer than 2 fragments contribute.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static float ComputeCosineAtTimePoint(
            ReadOnlySpan<float> matrix,
            ReadOnlySpan<float> libraryIntensities,
            int fragmentCount,
            int timePoint)
        {
            float dot = 0f, normLib = 0f, normObs = 0f;
            int nContrib = 0;
            int rowOffset = timePoint * fragmentCount;

            int limit = Math.Min(fragmentCount, libraryIntensities.Length);
            for (int f = 0; f < limit; f++)
            {
                float obs = matrix[rowOffset + f];
                float lib = libraryIntensities[f];
                if (obs <= 0f || lib <= 0f) continue;

                dot += obs * lib;
                normObs += obs * obs;
                normLib += lib * lib;
                nContrib++;
            }

            if (nContrib < 2) return float.NaN;

            float denom = MathF.Sqrt(normLib) * MathF.Sqrt(normObs);
            if (denom <= 0f) return float.NaN;

            return Math.Clamp(dot / denom, -1f, 1f);
        }

        /// <summary>
        /// Counts nonzero intensity values in a specific fragment column across all time points.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static int CountNonzeroInColumn(
            ReadOnlySpan<float> matrix,
            int fragmentIndex,
            int fragmentCount,
            int timePointCount)
        {
            int count = 0;
            for (int t = 0; t < timePointCount; t++)
            {
                if (matrix[t * fragmentCount + fragmentIndex] > 0f)
                    count++;
            }
            return count;
        }
    }
}