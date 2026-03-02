// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Buffers;
using System.Runtime.CompilerServices;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Phase 13: Computes advanced discriminative features from the time × fragment matrix.
    /// 
    /// These features close the gap between our ~17-feature classifier and DIA-NN's ~73 features.
    /// All computations use data already in memory (the matrix built during assembly) — 
    /// no additional file reads required.
    /// 
    /// Feature categories implemented:
    ///   1. Best-Fragment Reference Curve (DIA-NN's core innovation for interference handling)
    ///   2. Per-Fragment Signal Ratio Deviation (interference detection without cross-precursor analysis)
    ///   3. Smoothed Fragment Correlations (noise reduction via min-of-3 smoothing)
    ///   4. Signal-to-Noise Ratio (explicit S/N calculation)
    ///   5. Log Total Intensity (raw signal strength)
    ///   6. Intensity Profile Scores (boundary signal ratios for peak shape quality)
    /// 
    /// Performance contract:
    ///   - No heap allocations in the hot path (uses stackalloc and ArrayPool)
    ///   - All O(F² × T) where F = fragments (~6-12), T = time points (~20-100)
    ///   - Designed to add &lt;1 second to the assembly step for 39K precursors
    /// </summary>
    public static class DiaFeatureCalculator
    {
        /// <summary>
        /// Computes all Phase 13 features for a single precursor from its time × fragment matrix.
        /// Call this after matrix construction and peak detection in AssembleResultsWithTemporalScoring.
        /// 
        /// The matrix layout is: matrix[t * fragmentCount + f] = intensity of fragment f at time t.
        /// </summary>
        /// <param name="matrix">Row-major time × fragment matrix (contiguous float span).</param>
        /// <param name="libraryIntensities">Library/predicted intensities per fragment (length = fragmentCount).</param>
        /// <param name="fragmentCount">Number of fragments (columns in matrix).</param>
        /// <param name="timePointCount">Number of time points (rows in matrix).</param>
        /// <param name="apexIndex">Index of the apex time point (from peak detection or full-window max).</param>
        /// <param name="peakLeftIndex">Left boundary of detected peak (-1 if no peak detected, uses full window).</param>
        /// <param name="peakRightIndex">Right boundary of detected peak (-1 if no peak detected, uses full window).</param>
        /// <param name="result">DiaSearchResult to populate with computed features.</param>
        public static void ComputeAllFeatures(
            ReadOnlySpan<float> matrix,
            ReadOnlySpan<float> libraryIntensities,
            int fragmentCount,
            int timePointCount,
            int apexIndex,
            int peakLeftIndex,
            int peakRightIndex,
            DiaSearchResult result)
        {
            if (fragmentCount < 2 || timePointCount < 3)
            {
                SetNaNDefaults(result);
                return;
            }

            // Use peak boundaries if available, otherwise full window
            int rangeStart = peakLeftIndex >= 0 ? peakLeftIndex : 0;
            int rangeEnd = peakRightIndex >= 0 ? peakRightIndex : timePointCount - 1;
            int rangeLength = rangeEnd - rangeStart + 1;

            if (rangeLength < 3)
            {
                SetNaNDefaults(result);
                return;
            }

            // ═══════════════════════════════════════════════════════════════════
            //  FEATURE 1: Best-Fragment Reference Curve
            //  
            //  DIA-NN's core innovation: find the fragment whose XIC has the
            //  highest sum of Pearson correlations with all other fragments.
            //  This fragment is least affected by interference and serves as
            //  the "gold standard" elution profile.
            // ═══════════════════════════════════════════════════════════════════
            ComputeBestFragmentFeatures(matrix, libraryIntensities, fragmentCount,
                timePointCount, rangeStart, rangeEnd, apexIndex, result);

            // ═══════════════════════════════════════════════════════════════════
            //  FEATURE 2: Per-Fragment Signal Ratio Deviation
            //  
            //  Compares observed vs library intensity ratios per fragment at apex.
            //  True peptides match library ratios; interference inflates specific
            //  fragments, creating ratio deviations.
            // ═══════════════════════════════════════════════════════════════════
            ComputeSignalRatioFeatures(matrix, libraryIntensities, fragmentCount,
                apexIndex, result);

            // ═══════════════════════════════════════════════════════════════════
            //  FEATURE 3: Smoothed Fragment Correlations
            //  
            //  DIA-NN computes correlations on min-of-3-consecutive smoothed XICs.
            //  This reduces the impact of single-scan noise spikes that can inflate
            //  or deflate pairwise correlations.
            // ═══════════════════════════════════════════════════════════════════
            ComputeSmoothedCorrelationFeatures(matrix, fragmentCount, timePointCount,
                rangeStart, rangeEnd, result);

            // ═══════════════════════════════════════════════════════════════════
            //  FEATURE 4: Signal-to-Noise Ratio
            // ═══════════════════════════════════════════════════════════════════
            ComputeSignalToNoise(matrix, fragmentCount, timePointCount,
                rangeStart, rangeEnd, apexIndex, result);

            // ═══════════════════════════════════════════════════════════════════
            //  FEATURE 5: Log Total Intensity
            // ═══════════════════════════════════════════════════════════════════
            float totalIntensity = 0f;
            for (int f = 0; f < fragmentCount; f++)
            {
                for (int t = rangeStart; t <= rangeEnd; t++)
                {
                    totalIntensity += matrix[t * fragmentCount + f];
                }
            }
            result.LogTotalIntensity = totalIntensity > 0f ? MathF.Log2(totalIntensity) : 0f;

            // ═══════════════════════════════════════════════════════════════════
            //  FEATURE 6: Intensity Profile / Peak Shape Quality
            // ═══════════════════════════════════════════════════════════════════
            ComputePeakShapeFeatures(matrix, fragmentCount, timePointCount,
                rangeStart, rangeEnd, apexIndex, result);
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Best-Fragment Reference Curve
        // ─────────────────────────────────────────────────────────────────────

        private static void ComputeBestFragmentFeatures(
            ReadOnlySpan<float> matrix,
            ReadOnlySpan<float> libraryIntensities,
            int fragmentCount,
            int timePointCount,
            int rangeStart,
            int rangeEnd,
            int apexIndex,
            DiaSearchResult result)
        {
            int rangeLength = rangeEnd - rangeStart + 1;

            // Identify fragments with sufficient data points in the range
            Span<int> activeFrags = stackalloc int[Math.Min(fragmentCount, 64)];
            int nActive = 0;
            for (int f = 0; f < fragmentCount && nActive < activeFrags.Length; f++)
            {
                int nonzero = 0;
                for (int t = rangeStart; t <= rangeEnd; t++)
                {
                    if (matrix[t * fragmentCount + f] > 0f) nonzero++;
                }
                if (nonzero >= 3) activeFrags[nActive++] = f;
            }

            if (nActive < 2)
            {
                result.BestFragCorrelationSum = float.NaN;
                result.MedianFragRefCorr = float.NaN;
                result.MinFragRefCorr = float.NaN;
                result.StdFragRefCorr = float.NaN;
                return;
            }

            // Compute pairwise correlation matrix (upper triangle)
            // corrMatrix[a * nActive + b] stores correlation between activeFrags[a] and activeFrags[b]
            int corrSize = nActive * nActive;
            float[] corrRented = ArrayPool<float>.Shared.Rent(corrSize);
            Span<float> corrMatrix = corrRented.AsSpan(0, corrSize);

            try
            {
                // Fill correlation matrix
                for (int a = 0; a < nActive; a++)
                {
                    corrMatrix[a * nActive + a] = 1.0f; // self-correlation
                    for (int b = a + 1; b < nActive; b++)
                    {
                        float r = PearsonOnRange(matrix, activeFrags[a], activeFrags[b],
                            fragmentCount, rangeStart, rangeEnd);
                        if (float.IsNaN(r)) r = 0f;
                        corrMatrix[a * nActive + b] = r;
                        corrMatrix[b * nActive + a] = r;
                    }
                }

                // Find best fragment: max sum of correlations with all others
                int bestIdx = 0;
                float bestSum = float.MinValue;
                for (int a = 0; a < nActive; a++)
                {
                    float sum = 0f;
                    for (int b = 0; b < nActive; b++)
                    {
                        if (a != b) sum += corrMatrix[a * nActive + b];
                    }
                    if (sum > bestSum) { bestSum = sum; bestIdx = a; }
                }

                result.BestFragCorrelationSum = bestSum;
                result.BestFragIndex = activeFrags[bestIdx];

                // Compute individual fragment correlations with the best fragment
                Span<float> refCorrs = stackalloc float[nActive];
                int nRefCorrs = 0;
                for (int a = 0; a < nActive; a++)
                {
                    if (a == bestIdx) continue;
                    refCorrs[nRefCorrs++] = corrMatrix[a * nActive + bestIdx];
                }

                if (nRefCorrs > 0)
                {
                    // Sort for median
                    Span<float> sortable = refCorrs.Slice(0, nRefCorrs);
                    sortable.Sort();

                    result.MinFragRefCorr = sortable[0];
                    result.MedianFragRefCorr = nRefCorrs % 2 == 0
                        ? (sortable[nRefCorrs / 2 - 1] + sortable[nRefCorrs / 2]) / 2f
                        : sortable[nRefCorrs / 2];

                    // Standard deviation
                    float mean = 0f;
                    for (int i = 0; i < nRefCorrs; i++) mean += sortable[i];
                    mean /= nRefCorrs;
                    float variance = 0f;
                    for (int i = 0; i < nRefCorrs; i++)
                    {
                        float diff = sortable[i] - mean;
                        variance += diff * diff;
                    }
                    result.StdFragRefCorr = MathF.Sqrt(variance / nRefCorrs);
                }

                // Best-fragment weighted cosine at apex
                // Weight each fragment's contribution by its correlation with the best fragment
                if (apexIndex >= rangeStart && apexIndex <= rangeEnd)
                {
                    float dot = 0f, normLib = 0f, normObs = 0f;
                    int rowOffset = apexIndex * fragmentCount;
                    for (int a = 0; a < nActive; a++)
                    {
                        int f = activeFrags[a];
                        float obs = matrix[rowOffset + f];
                        if (obs <= 0f) continue;
                        float lib = f < libraryIntensities.Length ? libraryIntensities[f] : 0f;
                        float weight = a == bestIdx ? 1.0f : Math.Max(0f, corrMatrix[a * nActive + bestIdx]);
                        float wObs = obs * weight;
                        float wLib = lib * weight;
                        dot += wLib * wObs;
                        normLib += wLib * wLib;
                        normObs += wObs * wObs;
                    }
                    result.BestFragWeightedCosine = (normLib > 0f && normObs > 0f)
                        ? Math.Clamp(dot / (MathF.Sqrt(normLib) * MathF.Sqrt(normObs)), 0f, 1f)
                        : float.NaN;
                }
            }
            finally
            {
                ArrayPool<float>.Shared.Return(corrRented);
            }
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Per-Fragment Signal Ratio Deviation
        // ─────────────────────────────────────────────────────────────────────

        private static void ComputeSignalRatioFeatures(
            ReadOnlySpan<float> matrix,
            ReadOnlySpan<float> libraryIntensities,
            int fragmentCount,
            int apexIndex,
            DiaSearchResult result)
        {
            int rowOffset = apexIndex * fragmentCount;

            // Compute observed total and library total at apex
            float obsTotal = 0f;
            float libTotal = 0f;
            int activeCount = 0;
            for (int f = 0; f < fragmentCount; f++)
            {
                float obs = matrix[rowOffset + f];
                float lib = f < libraryIntensities.Length ? libraryIntensities[f] : 0f;
                if (obs > 0f && lib > 0f)
                {
                    obsTotal += obs;
                    libTotal += lib;
                    activeCount++;
                }
            }

            if (activeCount < 2 || obsTotal <= 0f || libTotal <= 0f)
            {
                result.MeanSignalRatioDeviation = float.NaN;
                result.MaxSignalRatioDeviation = float.NaN;
                result.StdSignalRatioDeviation = float.NaN;
                return;
            }

            // For each fragment: compute log-ratio of (obs/obsTotal) vs (lib/libTotal)
            // True peptides: ratios match → log-ratio ≈ 0
            // Interference: specific fragments inflated → large |log-ratio|
            Span<float> deviations = stackalloc float[Math.Min(fragmentCount, 64)];
            int nDev = 0;

            for (int f = 0; f < fragmentCount && nDev < deviations.Length; f++)
            {
                float obs = matrix[rowOffset + f];
                float lib = f < libraryIntensities.Length ? libraryIntensities[f] : 0f;
                if (obs <= 0f || lib <= 0f) continue;

                float obsRatio = obs / obsTotal;
                float libRatio = lib / libTotal;
                float logRatio = MathF.Abs(MathF.Log2(obsRatio / libRatio));
                deviations[nDev++] = logRatio;
            }

            if (nDev < 2)
            {
                result.MeanSignalRatioDeviation = float.NaN;
                result.MaxSignalRatioDeviation = float.NaN;
                result.StdSignalRatioDeviation = float.NaN;
                return;
            }

            float sum = 0f, maxDev = 0f;
            for (int i = 0; i < nDev; i++)
            {
                sum += deviations[i];
                if (deviations[i] > maxDev) maxDev = deviations[i];
            }
            float mean = sum / nDev;

            float varSum = 0f;
            for (int i = 0; i < nDev; i++)
            {
                float diff = deviations[i] - mean;
                varSum += diff * diff;
            }

            result.MeanSignalRatioDeviation = mean;
            result.MaxSignalRatioDeviation = maxDev;
            result.StdSignalRatioDeviation = MathF.Sqrt(varSum / nDev);
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Smoothed Fragment Correlations (min-of-3 consecutive)
        // ─────────────────────────────────────────────────────────────────────

        private static void ComputeSmoothedCorrelationFeatures(
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int timePointCount,
            int rangeStart,
            int rangeEnd,
            DiaSearchResult result)
        {
            int rangeLength = rangeEnd - rangeStart + 1;
            if (rangeLength < 5) // Need at least 5 time points for meaningful smoothing
            {
                result.SmoothedMeanFragCorr = float.NaN;
                result.SmoothedMinFragCorr = float.NaN;
                return;
            }

            // Build smoothed matrix using min-of-3 consecutive values per fragment
            int smoothedLength = rangeLength - 2; // lose one point on each end
            int smoothedSize = smoothedLength * fragmentCount;
            float[] smoothedRented = ArrayPool<float>.Shared.Rent(smoothedSize);
            Span<float> smoothed = smoothedRented.AsSpan(0, smoothedSize);

            try
            {
                for (int f = 0; f < fragmentCount; f++)
                {
                    for (int s = 0; s < smoothedLength; s++)
                    {
                        int t = rangeStart + s + 1; // center of the 3-point window
                        float v0 = matrix[(t - 1) * fragmentCount + f];
                        float v1 = matrix[t * fragmentCount + f];
                        float v2 = matrix[(t + 1) * fragmentCount + f];
                        smoothed[s * fragmentCount + f] = Math.Min(v0, Math.Min(v1, v2));
                    }
                }

                // Find active fragments in smoothed data
                Span<int> activeFrags = stackalloc int[Math.Min(fragmentCount, 64)];
                int nActive = 0;
                for (int f = 0; f < fragmentCount && nActive < activeFrags.Length; f++)
                {
                    int nonzero = 0;
                    for (int s = 0; s < smoothedLength; s++)
                    {
                        if (smoothed[s * fragmentCount + f] > 0f) nonzero++;
                    }
                    if (nonzero >= 3) activeFrags[nActive++] = f;
                }

                if (nActive < 2)
                {
                    result.SmoothedMeanFragCorr = float.NaN;
                    result.SmoothedMinFragCorr = float.NaN;
                    return;
                }

                // Pairwise correlations on smoothed data
                float sumCorr = 0f;
                float minCorr = float.MaxValue;
                int nPairs = 0;

                for (int a = 0; a < nActive; a++)
                {
                    for (int b = a + 1; b < nActive; b++)
                    {
                        float r = PearsonOnSmoothed(smoothed, activeFrags[a], activeFrags[b],
                            fragmentCount, smoothedLength);
                        if (!float.IsNaN(r))
                        {
                            sumCorr += r;
                            if (r < minCorr) minCorr = r;
                            nPairs++;
                        }
                    }
                }

                result.SmoothedMeanFragCorr = nPairs > 0 ? sumCorr / nPairs : float.NaN;
                result.SmoothedMinFragCorr = nPairs > 0 ? minCorr : float.NaN;
            }
            finally
            {
                ArrayPool<float>.Shared.Return(smoothedRented);
            }
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Signal-to-Noise Ratio
        // ─────────────────────────────────────────────────────────────────────

        private static void ComputeSignalToNoise(
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int timePointCount,
            int rangeStart,
            int rangeEnd,
            int apexIndex,
            DiaSearchResult result)
        {
            // Signal: total intensity at apex scan
            float signalAtApex = 0f;
            int apexRow = apexIndex * fragmentCount;
            for (int f = 0; f < fragmentCount; f++)
                signalAtApex += matrix[apexRow + f];

            if (signalAtApex <= 0f)
            {
                result.Log2SignalToNoise = float.NaN;
                return;
            }

            // Noise: median of per-time-point total intensities OUTSIDE the peak range
            // If no peak detected (rangeStart=0, rangeEnd=last), use the bottom quartile
            Span<float> noiseValues = stackalloc float[Math.Min(timePointCount, 256)];
            int nNoise = 0;

            for (int t = 0; t < timePointCount && nNoise < noiseValues.Length; t++)
            {
                // Use points outside the peak if we have peak boundaries
                if (t >= rangeStart && t <= rangeEnd) continue;

                float total = 0f;
                int rowOff = t * fragmentCount;
                for (int f = 0; f < fragmentCount; f++)
                    total += matrix[rowOff + f];
                noiseValues[nNoise++] = total;
            }

            if (nNoise < 2)
            {
                // Not enough noise points — use a rough estimate:
                // take the 25th percentile of all time points
                nNoise = 0;
                for (int t = 0; t < timePointCount && nNoise < noiseValues.Length; t++)
                {
                    float total = 0f;
                    int rowOff = t * fragmentCount;
                    for (int f = 0; f < fragmentCount; f++)
                        total += matrix[rowOff + f];
                    noiseValues[nNoise++] = total;
                }
                if (nNoise == 0)
                {
                    result.Log2SignalToNoise = float.NaN;
                    return;
                }
                Span<float> allSorted = noiseValues.Slice(0, nNoise);
                allSorted.Sort();
                int q25Idx = nNoise / 4;
                float noiseEstimate = allSorted[q25Idx];
                result.Log2SignalToNoise = noiseEstimate > 0f
                    ? MathF.Log2(signalAtApex / noiseEstimate)
                    : MathF.Log2(signalAtApex); // If noise is zero, just use signal magnitude
                return;
            }

            Span<float> sorted = noiseValues.Slice(0, nNoise);
            sorted.Sort();
            float medianNoise = nNoise % 2 == 0
                ? (sorted[nNoise / 2 - 1] + sorted[nNoise / 2]) / 2f
                : sorted[nNoise / 2];

            result.Log2SignalToNoise = medianNoise > 0f
                ? MathF.Log2(signalAtApex / medianNoise)
                : MathF.Log2(signalAtApex);
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Peak Shape Features
        // ─────────────────────────────────────────────────────────────────────

        private static void ComputePeakShapeFeatures(
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int timePointCount,
            int rangeStart,
            int rangeEnd,
            int apexIndex,
            DiaSearchResult result)
        {
            // Signal at apex
            float apexSignal = 0f;
            int apexRow = apexIndex * fragmentCount;
            for (int f = 0; f < fragmentCount; f++)
                apexSignal += matrix[apexRow + f];

            if (apexSignal <= 0f)
            {
                result.BoundarySignalRatio = float.NaN;
                result.ApexToMeanRatio = float.NaN;
                return;
            }

            // Signal at boundaries (average of left and right boundary)
            float leftSignal = 0f, rightSignal = 0f;
            int leftRow = rangeStart * fragmentCount;
            int rightRow = rangeEnd * fragmentCount;
            for (int f = 0; f < fragmentCount; f++)
            {
                leftSignal += matrix[leftRow + f];
                rightSignal += matrix[rightRow + f];
            }
            float boundaryAvg = (leftSignal + rightSignal) / 2f;

            // BoundarySignalRatio: low for true peaks (boundaries should be near baseline),
            // high for interference (signal doesn't drop at boundaries)
            result.BoundarySignalRatio = boundaryAvg / apexSignal;

            // ApexToMeanRatio: how much the apex stands above the average
            float totalSignal = 0f;
            int rangeLength = rangeEnd - rangeStart + 1;
            for (int t = rangeStart; t <= rangeEnd; t++)
            {
                int rowOff = t * fragmentCount;
                for (int f = 0; f < fragmentCount; f++)
                    totalSignal += matrix[rowOff + f];
            }
            float meanSignal = totalSignal / rangeLength;
            result.ApexToMeanRatio = meanSignal > 0f ? apexSignal / meanSignal : 1f;
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Helper: Pearson correlation between two fragments over a range
        // ─────────────────────────────────────────────────────────────────────

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static float PearsonOnRange(
            ReadOnlySpan<float> matrix, int fragA, int fragB,
            int fragmentCount, int rangeStart, int rangeEnd)
        {
            float sumA = 0, sumB = 0, sumAB = 0, sumA2 = 0, sumB2 = 0;
            int n = 0;

            for (int t = rangeStart; t <= rangeEnd; t++)
            {
                float a = matrix[t * fragmentCount + fragA];
                float b = matrix[t * fragmentCount + fragB];
                if (a <= 0f || b <= 0f) continue;

                sumA += a; sumB += b;
                sumAB += a * b;
                sumA2 += a * a; sumB2 += b * b;
                n++;
            }

            if (n < 3) return float.NaN;

            float denom = (n * sumA2 - sumA * sumA) * (n * sumB2 - sumB * sumB);
            if (denom <= 0f) return float.NaN;

            return Math.Clamp((n * sumAB - sumA * sumB) / MathF.Sqrt(denom), -1f, 1f);
        }

        /// <summary>
        /// Pearson correlation on the smoothed matrix (different layout — all points used).
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static float PearsonOnSmoothed(
            ReadOnlySpan<float> smoothed, int fragA, int fragB,
            int fragmentCount, int smoothedLength)
        {
            float sumA = 0, sumB = 0, sumAB = 0, sumA2 = 0, sumB2 = 0;
            int n = 0;

            for (int s = 0; s < smoothedLength; s++)
            {
                float a = smoothed[s * fragmentCount + fragA];
                float b = smoothed[s * fragmentCount + fragB];
                if (a <= 0f || b <= 0f) continue;

                sumA += a; sumB += b;
                sumAB += a * b;
                sumA2 += a * a; sumB2 += b * b;
                n++;
            }

            if (n < 3) return float.NaN;

            float denom = (n * sumA2 - sumA * sumA) * (n * sumB2 - sumB * sumB);
            if (denom <= 0f) return float.NaN;

            return Math.Clamp((n * sumAB - sumA * sumB) / MathF.Sqrt(denom), -1f, 1f);
        }

        private static void SetNaNDefaults(DiaSearchResult result)
        {
            result.BestFragCorrelationSum = float.NaN;
            result.MedianFragRefCorr = float.NaN;
            result.MinFragRefCorr = float.NaN;
            result.StdFragRefCorr = float.NaN;
            result.BestFragWeightedCosine = float.NaN;
            result.BestFragIndex = -1;
            result.MeanSignalRatioDeviation = float.NaN;
            result.MaxSignalRatioDeviation = float.NaN;
            result.StdSignalRatioDeviation = float.NaN;
            result.SmoothedMeanFragCorr = float.NaN;
            result.SmoothedMinFragCorr = float.NaN;
            result.Log2SignalToNoise = float.NaN;
            result.LogTotalIntensity = 0f;
            result.BoundarySignalRatio = float.NaN;
            result.ApexToMeanRatio = float.NaN;
        }
    }
}