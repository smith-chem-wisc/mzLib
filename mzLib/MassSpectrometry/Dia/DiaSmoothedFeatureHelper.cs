// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Buffers;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Phase 13, Prompt 6: Smoothed correlation features and signal-to-noise.
    /// 
    /// ComputeSmoothedCorrelationFeatures:
    ///   Applies min-of-3 temporal smoothing to each fragment XIC (on a copy),
    ///   then computes pairwise Pearson correlations among detected fragments.
    ///   This noise-robust version of MeanFragCorr/MinFragCorr reduces the impact
    ///   of single-scan spikes or dropouts.
    /// 
    /// ComputeSignalToNoise:
    ///   Computes log2 signal-to-noise ratio by comparing total fragment intensity
    ///   at the apex scan versus the median total intensity across all scans.
    /// </summary>
    internal static class DiaSmoothedFeatureHelper
    {
        /// <summary>
        /// Computes smoothed pairwise fragment correlations and populates the result.
        /// 
        /// Algorithm:
        ///   1. Copy the intensity matrix (to avoid modifying the original)
        ///   2. Apply min-of-3 temporal smoothing to each fragment XIC:
        ///      smoothed[t] = min(original[t-1], original[t], original[t+1])
        ///      This removes single-scan spikes while preserving real signal
        ///   3. Find detected fragments (≥3 nonzero time points after smoothing)
        ///   4. Compute pairwise Pearson correlations on the smoothed data
        ///   5. Report mean and min correlation
        /// 
        /// Requires ≥2 detected fragments with ≥3 time points each.
        /// </summary>
        /// <param name="matrix">
        /// Flat row-major intensity matrix [timePointCount × fragmentCount].
        /// NOT modified — the method works on an internal copy.
        /// </param>
        /// <param name="fragmentCount">Number of fragment ions (columns).</param>
        /// <param name="timePointCount">Number of time points (rows).</param>
        /// <param name="result">DiaSearchResult to populate with SmoothedMeanFragCorr, SmoothedMinFragCorr.</param>
        public static void ComputeSmoothedCorrelationFeatures(
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int timePointCount,
            DiaSearchResult result)
        {
            if (fragmentCount < 2 || timePointCount < 3)
                return;

            int matrixSize = timePointCount * fragmentCount;

            // Copy the matrix for smoothing (pooled to avoid per-precursor GC pressure)
            float[] smoothedRented = ArrayPool<float>.Shared.Rent(matrixSize);
            try
            {
                Array.Clear(smoothedRented, 0, matrixSize);
                matrix.Slice(0, matrixSize).CopyTo(smoothedRented.AsSpan(0, matrixSize));
                float[] smoothed = smoothedRented;

                // Apply min-of-3 temporal smoothing per fragment column
                // smoothed[t,f] = min(original[t-1,f], original[t,f], original[t+1,f])
                // Process in-place using a rolling approach per column
                for (int f = 0; f < fragmentCount; f++)
                {
                    // We need the original values, so process with a temp buffer per column
                    float prev = matrix[0 * fragmentCount + f];
                    float curr = matrix[0 * fragmentCount + f];

                    for (int t = 0; t < timePointCount; t++)
                    {
                        float next = (t + 1 < timePointCount)
                            ? matrix[(t + 1) * fragmentCount + f]
                            : matrix[t * fragmentCount + f];

                        // At boundaries: t=0 uses min(curr, next), t=last uses min(prev, curr)
                        float minVal;
                        if (t == 0)
                            minVal = MathF.Min(curr, next);
                        else if (t == timePointCount - 1)
                            minVal = MathF.Min(prev, curr);
                        else
                            minVal = MathF.Min(prev, MathF.Min(curr, next));

                        smoothed[t * fragmentCount + f] = minVal;

                        prev = curr;
                        curr = next;
                    }
                }

                // Find detected fragments (≥3 nonzero time points after smoothing)
                Span<int> detectedIndices = stackalloc int[Math.Min(fragmentCount, 64)];
                int nDetected = 0;

                for (int f = 0; f < fragmentCount && nDetected < detectedIndices.Length; f++)
                {
                    int nonzero = 0;
                    for (int t = 0; t < timePointCount; t++)
                    {
                        if (smoothed[t * fragmentCount + f] > 0f)
                            nonzero++;
                    }
                    if (nonzero >= 3)
                        detectedIndices[nDetected++] = f;
                }

                if (nDetected < 2)
                    return;

                // Compute pairwise Pearson correlations on smoothed data
                float sumCorr = 0f;
                float minCorr = float.MaxValue;
                int nPairs = 0;

                ReadOnlySpan<float> smoothedSpan = smoothed.AsSpan();

                for (int a = 0; a < nDetected; a++)
                {
                    for (int b = a + 1; b < nDetected; b++)
                    {
                        float r = PearsonCorrelation(
                            smoothedSpan, detectedIndices[a], detectedIndices[b],
                            fragmentCount, timePointCount);

                        if (!float.IsNaN(r))
                        {
                            sumCorr += r;
                            if (r < minCorr) minCorr = r;
                            nPairs++;
                        }
                    }
                }

                if (nPairs > 0)
                {
                    result.SmoothedMeanFragCorr = sumCorr / nPairs;
                    result.SmoothedMinFragCorr = minCorr;
                }
            }
            finally
            {
                ArrayPool<float>.Shared.Return(smoothedRented);
            }
        }

        /// <summary>
        /// Computes log2 signal-to-noise ratio and populates the result.
        /// 
        /// Algorithm:
        ///   1. For each time point, compute total fragment intensity (TIC)
        ///   2. Signal = TIC at the apex scan
        ///   3. Noise = median TIC across all scans (robust background estimate)
        ///   4. S/N = signal / noise (or signal / small_value if noise ≈ 0)
        ///   5. Report log2(S/N)
        /// 
        /// Requires ≥3 time points. Returns NaN if apex signal is zero.
        /// </summary>
        /// <param name="matrix">
        /// Flat row-major intensity matrix [timePointCount × fragmentCount].
        /// </param>
        /// <param name="apexScanIndex">Local apex scan index (row in the matrix).</param>
        /// <param name="fragmentCount">Number of fragment ions (columns).</param>
        /// <param name="timePointCount">Number of time points (rows).</param>
        /// <param name="result">DiaSearchResult to populate with Log2SignalToNoise.</param>
        public static void ComputeSignalToNoise(
            ReadOnlySpan<float> matrix,
            int apexScanIndex,
            int fragmentCount,
            int timePointCount,
            DiaSearchResult result)
        {
            if (timePointCount < 3 || apexScanIndex < 0 || apexScanIndex >= timePointCount)
                return;

            // Compute TIC per time point
            Span<float> tics = timePointCount <= 256
                ? stackalloc float[timePointCount]
                : new float[timePointCount];

            for (int t = 0; t < timePointCount; t++)
            {
                float total = 0f;
                int rowOffset = t * fragmentCount;
                for (int f = 0; f < fragmentCount; f++)
                    total += matrix[rowOffset + f];
                tics[t] = total;
            }

            float apexTic = tics[apexScanIndex];
            if (apexTic <= 0f)
                return;

            // Find median TIC (sort a copy)
            Span<float> sorted = timePointCount <= 256
                ? stackalloc float[timePointCount]
                : new float[timePointCount];
            tics.CopyTo(sorted);
            SortSpan(sorted);

            float medianTic = timePointCount % 2 == 1
                ? sorted[timePointCount / 2]
                : (sorted[timePointCount / 2 - 1] + sorted[timePointCount / 2]) / 2f;

            // Avoid division by zero: use a small floor for noise
            const float noiseFloor = 1f;
            float noise = MathF.Max(medianTic, noiseFloor);

            result.Log2SignalToNoise = MathF.Log2(apexTic / noise);
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
        /// Simple insertion sort for small spans (typically ≤ 30 elements).
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