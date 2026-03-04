// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Phase 13, Prompt 4: Best-fragment reference curve features.
    /// Phase 16A, Prompt 1: BestFragWeightedCosine confirmed active in classifier (feature [26]).
    /// 
    /// Identifies the "best" (least interfered) fragment ion by finding the one whose
    /// XIC has the highest sum of Pearson correlations with all other detected fragments.
    /// Then computes each detected fragment's correlation with this reference fragment,
    /// producing summary statistics (sum, median, min, std) that measure how consistently
    /// all fragments coelute with the cleanest signal.
    /// 
    /// Phase 13 Action Item 1: Also computes BestFragWeightedCosine — a cosine score
    /// between library and observed apex intensities, weighted by each fragment's
    /// correlation with the best fragment. This down-weights interfered fragments.
    /// 
    /// True peptide signals show high, consistent correlations with the best fragment.
    /// Interfered or random matches show low or scattered correlations.
    /// </summary>
    internal static class DiaBestFragmentHelper
    {
        /// <summary>
        /// Computes best-fragment reference curve features and populates the result.
        /// 
        /// Algorithm:
        ///   1. Identify detected fragments (≥3 nonzero time points)
        ///   2. Compute full pairwise Pearson correlation matrix among detected fragments
        ///   3. Find the "best" fragment = the one with the highest sum of correlations
        ///   4. Extract each fragment's correlation with the best fragment
        ///   5. Compute summary stats: sum, median, min, std of those correlations
        ///   6. Compute weighted cosine at apex using correlations as weights (Action Item 1)
        /// 
        /// Requires ≥3 detected fragments to produce meaningful results.
        /// </summary>
        /// <param name="matrix">
        /// Flat row-major intensity matrix [timePointCount × fragmentCount].
        /// Access pattern: matrix[t * fragmentCount + f].
        /// </param>
        /// <param name="fragmentCount">Number of fragment ions (columns).</param>
        /// <param name="timePointCount">Number of time points (rows).</param>
        /// <param name="detectedFragmentCount">
        /// Number of fragments with at least one nonzero data point.
        /// Used as an early-exit check.
        /// </param>
        /// <param name="result">DiaSearchResult to populate with computed features.</param>
        /// <param name="libraryIntensities">
        /// Library-predicted fragment intensities (parallel to fragment columns).
        /// Used for weighted cosine computation.
        /// </param>
        /// <param name="apexRowOffset">
        /// Offset into the matrix for the apex scan row (= apexScanIndex * fragmentCount).
        /// Used for weighted cosine computation.
        /// </param>
        public static void ComputeBestFragmentFeatures(
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int timePointCount,
            int detectedFragmentCount,
            DiaSearchResult result,
            ReadOnlySpan<float> libraryIntensities,
            int apexRowOffset)
        {
            // Need at least 3 detected fragments for meaningful pairwise correlations
            if (detectedFragmentCount < 3 || fragmentCount < 3 || timePointCount < 3)
                return; // fields remain NaN from constructor

            // Step 1: Find fragments with ≥3 nonzero time points (suitable for Pearson)
            Span<int> detectedIndices = stackalloc int[Math.Min(fragmentCount, 64)];
            int nDetected = 0;

            for (int f = 0; f < fragmentCount && nDetected < detectedIndices.Length; f++)
            {
                int nonzero = 0;
                for (int t = 0; t < timePointCount; t++)
                {
                    if (matrix[t * fragmentCount + f] > 0f)
                        nonzero++;
                }
                if (nonzero >= 3)
                    detectedIndices[nDetected++] = f;
            }

            if (nDetected < 3)
                return;

            // Step 2: Compute pairwise Pearson correlation matrix (upper triangle)
            // Store as flat array: corrMatrix[i * nDetected + j]
            // We only need the upper triangle but store symmetrically for easy row-sum
            int corrSize = nDetected * nDetected;
            Span<float> corrMatrix = corrSize <= 256
                ? stackalloc float[corrSize]
                : new float[corrSize];
            corrMatrix.Clear();

            for (int a = 0; a < nDetected; a++)
            {
                corrMatrix[a * nDetected + a] = 1.0f; // self-correlation
                for (int b = a + 1; b < nDetected; b++)
                {
                    float r = PearsonCorrelation(
                        matrix, detectedIndices[a], detectedIndices[b],
                        fragmentCount, timePointCount);

                    if (float.IsNaN(r))
                        r = 0f; // treat undefined correlation as zero

                    corrMatrix[a * nDetected + b] = r;
                    corrMatrix[b * nDetected + a] = r;
                }
            }

            // Step 3: Find the best fragment (highest sum of correlations with others)
            int bestIdx = 0;
            float bestSum = float.MinValue;

            for (int a = 0; a < nDetected; a++)
            {
                float rowSum = 0f;
                for (int b = 0; b < nDetected; b++)
                {
                    if (b != a)
                        rowSum += corrMatrix[a * nDetected + b];
                }
                if (rowSum > bestSum)
                {
                    bestSum = rowSum;
                    bestIdx = a;
                }
            }

            // Step 4: Extract each fragment's correlation with the best fragment
            // (excluding the best fragment's self-correlation)
            Span<float> refCorrs = stackalloc float[nDetected - 1];
            int idx = 0;
            for (int a = 0; a < nDetected; a++)
            {
                if (a == bestIdx) continue;
                refCorrs[idx++] = corrMatrix[bestIdx * nDetected + a];
            }

            int n = refCorrs.Length;

            // Step 5: Compute summary statistics
            // BestFragCorrelationSum = sum of correlations between best and all others
            result.BestFragCorrelationSum = bestSum;

            // Sort for median
            SortSpan(refCorrs);

            // MedianFragRefCorr
            result.MedianFragRefCorr = n % 2 == 1
                ? refCorrs[n / 2]
                : (refCorrs[n / 2 - 1] + refCorrs[n / 2]) / 2f;

            // MinFragRefCorr (after sort, it's the first element)
            result.MinFragRefCorr = refCorrs[0];

            // StdFragRefCorr (population standard deviation)
            float sum = 0f;
            for (int i = 0; i < n; i++)
                sum += refCorrs[i];
            float mean = sum / n;

            float sumSqDev = 0f;
            for (int i = 0; i < n; i++)
            {
                float dev = refCorrs[i] - mean;
                sumSqDev += dev * dev;
            }
            result.StdFragRefCorr = MathF.Sqrt(sumSqDev / n);

            // Step 6: Best-fragment weighted cosine at apex (Action Item 1)
            // Weight each fragment's contribution to the cosine by its correlation
            // with the best fragment. Negative correlations are clamped to 0.
            // The best fragment itself gets weight 1.0.
            result.BestFragIndex = detectedIndices[bestIdx];
            ComputeWeightedCosine(
                matrix, corrMatrix, detectedIndices, nDetected, bestIdx,
                libraryIntensities, fragmentCount, apexRowOffset, result);
        }

        /// <summary>
        /// Computes weighted cosine between library and observed apex intensities,
        /// using each fragment's correlation with the best fragment as the weight.
        /// Fragments that correlate poorly with the reference are down-weighted,
        /// producing a cleaner cosine score that is robust to interference.
        /// </summary>
        private static void ComputeWeightedCosine(
            ReadOnlySpan<float> matrix,
            ReadOnlySpan<float> corrMatrix,
            ReadOnlySpan<int> detectedIndices,
            int nDetected,
            int bestIdx,
            ReadOnlySpan<float> libraryIntensities,
            int fragmentCount,
            int apexRowOffset,
            DiaSearchResult result)
        {
            float dot = 0f, normLib = 0f, normObs = 0f;

            for (int a = 0; a < nDetected; a++)
            {
                int f = detectedIndices[a];

                // Observed intensity at the apex scan for this fragment
                float obs = matrix[apexRowOffset + f];
                if (obs <= 0f) continue;

                // Library intensity for this fragment
                float lib = f < libraryIntensities.Length ? libraryIntensities[f] : 0f;
                if (lib <= 0f) continue;

                // Weight: best fragment gets 1.0, others get their correlation with best (clamped ≥ 0)
                float weight = (a == bestIdx)
                    ? 1.0f
                    : Math.Max(0f, corrMatrix[a * nDetected + bestIdx]);

                float wObs = obs * weight;
                float wLib = lib * weight;
                dot += wLib * wObs;
                normLib += wLib * wLib;
                normObs += wObs * wObs;
            }

            // Only set if we had valid weighted contributions
            if (normLib > 0f && normObs > 0f)
            {
                result.BestFragWeightedCosine = Math.Clamp(
                    dot / (MathF.Sqrt(normLib) * MathF.Sqrt(normObs)), 0f, 1f);
            }
            // else: remains NaN from constructor
        }

        /// <summary>
        /// Pearson correlation between two fragment XICs across all time points.
        /// Uses only time points where BOTH fragments have nonzero intensity.
        /// Returns NaN if fewer than 3 co-detected time points.
        /// </summary>
        private static float PearsonCorrelation(
            ReadOnlySpan<float> matrix, int fragA, int fragB,
            int fragmentCount, int timePointCount)
        {
            float sumA = 0, sumB = 0, sumAB = 0, sumA2 = 0, sumB2 = 0;
            int n = 0;

            for (int t = 0; t < timePointCount; t++)
            {
                float a = matrix[t * fragmentCount + fragA];
                float b = matrix[t * fragmentCount + fragB];
                if (a <= 0f || b <= 0f) continue;

                sumA += a;
                sumB += b;
                sumAB += a * b;
                sumA2 += a * a;
                sumB2 += b * b;
                n++;
            }

            if (n < 3) return float.NaN;

            float denom = (n * sumA2 - sumA * sumA) * (n * sumB2 - sumB * sumB);
            if (denom <= 0f) return float.NaN;

            float r = (n * sumAB - sumA * sumB) / MathF.Sqrt(denom);
            return Math.Clamp(r, -1f, 1f);
        }

        /// <summary>
        /// Simple insertion sort for small spans (typically ≤ 12 elements).
        /// </summary>
        private static void SortSpan(Span<float> span)
        {
            for (int i = 1; i < span.Length; i++)
            {
                float key = span[i];
                int j = i - 1;
                while (j >= 0 && span[j] > key)
                {
                    span[j + 1] = span[j];
                    j--;
                }
                span[j + 1] = key;
            }
        }
    }
}