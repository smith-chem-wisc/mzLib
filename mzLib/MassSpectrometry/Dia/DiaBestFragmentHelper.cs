// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Phase 13, Prompt 4: Best-fragment reference curve features.
    /// Phase 16A, Prompt 1: BestFragWeightedCosine confirmed active in classifier (feature [26]).
    /// Phase 19 (this revision): Pearson correlations now computed over an apex-centred
    /// window sized to the measured peak FWHM, replacing the full extraction window.
    ///
    /// Motivation:
    ///   The full extraction window (~50 scans) is dominated by baseline noise on either
    ///   side of the chromatographic peak. Correlating 50 points where ~45 are near-zero
    ///   produces spuriously high or noisy Pearson values. The apex-centred window
    ///   (typically 5-7 scans = ±FWHM/2 around the apex) captures only the region where
    ///   the peak actually lives, giving a physically meaningful correlation.
    ///
    /// Window sizing:
    ///   scanInterval   = median inter-scan gap in refRts
    ///   halfWidthScans = max(MinApexHalfWidthScans,
    ///                        round(peakWidthMinutes / 2 / scanInterval))
    ///   range          = [apexIdx - halfWidthScans, apexIdx + halfWidthScans]
    ///
    /// Fragment eligibility: ≥ 2 nonzero time points within the window
    ///   (relaxed from 3 because the window is narrow by design).
    ///
    /// Identifies the "best" (least interfered) fragment ion by finding the one whose
    /// XIC has the highest sum of Pearson correlations with all other detected fragments.
    /// Then computes each detected fragment's correlation with this reference fragment,
    /// producing summary statistics (sum, median, min, std) that measure how consistently
    /// all fragments coelute with the cleanest signal.
    ///
    /// Also computes BestFragWeightedCosine — a cosine score between library and observed
    /// apex intensities, weighted by each fragment's correlation with the best fragment.
    /// </summary>
    internal static class DiaBestFragmentHelper
    {
        /// <summary>
        /// Minimum half-width in scans on each side of the apex.
        /// Guarantees at least this many scans even for very narrow peaks.
        /// </summary>
        private const int MinApexHalfWidthScans = 3;

        /// <summary>
        /// Computes best-fragment reference curve features and populates the result.
        ///
        /// Algorithm:
        ///   1. Compute apex-centred window from peakWidthMinutes and refRts scan interval
        ///   2. Identify detected fragments (≥2 nonzero time points within window)
        ///   3. Compute full pairwise Pearson correlation matrix within window
        ///   4. Find the "best" fragment = highest sum of correlations with others
        ///   5. Extract each fragment's correlation with the best fragment
        ///   6. Compute summary stats: sum, median, min, std of those correlations
        ///   7. Compute weighted cosine at apex using correlations as weights
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
        /// Number of fragments with at least one nonzero data point across the full window.
        /// Used as an early-exit check only.
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
        /// <param name="apexLocalIdx">
        /// Index of the apex scan within the local refRts / matrix time axis.
        /// Centre of the correlation window.
        /// </param>
        /// <param name="peakWidthMinutes">
        /// Measured chromatographic peak width in minutes
        /// (refRts[RightIndex] - refRts[LeftIndex] from the detected peak group).
        /// Used to derive halfWidthScans. Pass 0f to use MinApexHalfWidthScans only.
        /// </param>
        /// <param name="refRts">
        /// RT timestamps of the reference fragment, length = timePointCount.
        /// Used to compute the median inter-scan interval for window sizing.
        /// </param>
        public static void ComputeBestFragmentFeatures(
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int timePointCount,
            int detectedFragmentCount,
            DiaSearchResult result,
            ReadOnlySpan<float> libraryIntensities,
            int apexRowOffset,
            int apexLocalIdx,
            float peakWidthMinutes,
            ReadOnlySpan<float> refRts)
        {
            // Need at least 3 detected fragments for meaningful pairwise correlations
            if (detectedFragmentCount < 3 || fragmentCount < 3 || timePointCount < 3)
                return;

            // ── Step 0: Compute apex-centred window ──────────────────────────
            int rangeStart, rangeEnd;
            ComputeApexWindow(timePointCount, apexLocalIdx, peakWidthMinutes, refRts,
                              out rangeStart, out rangeEnd);

            int rangeLength = rangeEnd - rangeStart + 1;
            if (rangeLength < 2) return;

            // ── Step 1: Find fragments with ≥2 nonzero points within window ──
            // Threshold relaxed from 3 to 2 because the window is narrow.
            Span<int> detectedIndices = stackalloc int[Math.Min(fragmentCount, 64)];
            int nDetected = 0;

            for (int f = 0; f < fragmentCount && nDetected < detectedIndices.Length; f++)
            {
                int nonzero = 0;
                for (int t = rangeStart; t <= rangeEnd; t++)
                    if (matrix[t * fragmentCount + f] > 0f) nonzero++;
                if (nonzero >= 2)
                    detectedIndices[nDetected++] = f;
            }

            if (nDetected < 3)
                return;

            // ── Step 2: Compute pairwise Pearson correlation matrix ───────────
            // Within the apex-centred window only.
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
                    float r = PearsonCorrelationOnRange(
                        matrix, detectedIndices[a], detectedIndices[b],
                        fragmentCount, rangeStart, rangeEnd);

                    if (float.IsNaN(r)) r = 0f; // undefined → zero

                    corrMatrix[a * nDetected + b] = r;
                    corrMatrix[b * nDetected + a] = r;
                }
            }

            // ── Step 3: Find best fragment (highest row-sum) ─────────────────
            int bestIdx = 0;
            float bestSum = float.MinValue;

            for (int a = 0; a < nDetected; a++)
            {
                float rowSum = 0f;
                for (int b = 0; b < nDetected; b++)
                    if (b != a) rowSum += corrMatrix[a * nDetected + b];
                if (rowSum > bestSum) { bestSum = rowSum; bestIdx = a; }
            }

            // ── Step 4: Extract correlations of all others vs best fragment ───
            Span<float> refCorrs = stackalloc float[nDetected - 1];
            int idx = 0;
            for (int a = 0; a < nDetected; a++)
            {
                if (a == bestIdx) continue;
                refCorrs[idx++] = corrMatrix[bestIdx * nDetected + a];
            }

            int n = refCorrs.Length;

            // ── Step 5: Summary statistics ────────────────────────────────────
            result.BestFragCorrelationSum = bestSum;

            SortSpan(refCorrs);

            result.MedianFragRefCorr = n % 2 == 1
                ? refCorrs[n / 2]
                : (refCorrs[n / 2 - 1] + refCorrs[n / 2]) / 2f;

            result.MinFragRefCorr = refCorrs[0];

            float sum = 0f;
            for (int i = 0; i < n; i++) sum += refCorrs[i];
            float mean = sum / n;
            float sumSqDev = 0f;
            for (int i = 0; i < n; i++) { float dev = refCorrs[i] - mean; sumSqDev += dev * dev; }
            result.StdFragRefCorr = MathF.Sqrt(sumSqDev / n);

            // ── Step 6: Best-fragment weighted cosine at apex ─────────────────
            result.BestFragIndex = detectedIndices[bestIdx];
            ComputeWeightedCosine(
                matrix, corrMatrix, detectedIndices, nDetected, bestIdx,
                libraryIntensities, fragmentCount, apexRowOffset, result);
        }

        // ── Window computation ────────────────────────────────────────────────

        /// <summary>
        /// Computes the apex-centred correlation window bounds.
        /// halfWidthScans = max(MinApexHalfWidthScans, round(peakWidthMinutes/2/scanInterval))
        /// </summary>
        private static void ComputeApexWindow(
            int timePointCount,
            int apexLocalIdx,
            float peakWidthMinutes,
            ReadOnlySpan<float> refRts,
            out int rangeStart,
            out int rangeEnd)
        {
            // Compute median inter-scan interval from refRts
            float scanInterval = 0.03f; // ~30ms fallback
            if (timePointCount >= 2)
            {
                // Collect gaps into a small buffer; use stack for short sequences
                if (timePointCount - 1 <= 64)
                {
                    Span<float> gaps = stackalloc float[timePointCount - 1];
                    for (int i = 0; i < timePointCount - 1; i++)
                        gaps[i] = refRts[i + 1] - refRts[i];
                    SortSpan(gaps);
                    float med = gaps[(timePointCount - 1) / 2];
                    if (med > 0f) scanInterval = med;
                }
                else
                {
                    // Fallback: global average for long sequences
                    float total = refRts[timePointCount - 1] - refRts[0];
                    float avg = total / (timePointCount - 1);
                    if (avg > 0f) scanInterval = avg;
                }
            }

            int halfWidthScans = MinApexHalfWidthScans;
            if (peakWidthMinutes > 0f)
            {
                int dynamic = (int)Math.Round(peakWidthMinutes / 2f / scanInterval);
                halfWidthScans = Math.Max(MinApexHalfWidthScans, dynamic);
            }

            rangeStart = Math.Max(0, apexLocalIdx - halfWidthScans);
            rangeEnd = Math.Min(timePointCount - 1, apexLocalIdx + halfWidthScans);
        }

        // ── Pearson correlation within a range ────────────────────────────────

        /// <summary>
        /// Pearson correlation between two fragment XICs within [rangeStart, rangeEnd].
        /// Uses only time points where BOTH fragments have nonzero intensity.
        /// Returns NaN if fewer than 2 co-detected time points.
        /// </summary>
        private static float PearsonCorrelationOnRange(
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
                sumA2 += a * a;
                sumB2 += b * b;
                n++;
            }

            if (n < 2) return float.NaN;

            float denom = (n * sumA2 - sumA * sumA) * (n * sumB2 - sumB * sumB);
            if (denom <= 0f) return float.NaN;

            float r = (n * sumAB - sumA * sumB) / MathF.Sqrt(denom);
            return Math.Clamp(r, -1f, 1f);
        }

        // ── Weighted cosine ───────────────────────────────────────────────────

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
                float obs = matrix[apexRowOffset + f];
                if (obs <= 0f) continue;
                float lib = f < libraryIntensities.Length ? libraryIntensities[f] : 0f;
                if (lib <= 0f) continue;

                float weight = (a == bestIdx)
                    ? 1.0f
                    : Math.Max(0f, corrMatrix[a * nDetected + bestIdx]);

                float wObs = obs * weight;
                float wLib = lib * weight;
                dot += wLib * wObs;
                normLib += wLib * wLib;
                normObs += wObs * wObs;
            }

            if (normLib > 0f && normObs > 0f)
                result.BestFragWeightedCosine = Math.Clamp(
                    dot / (MathF.Sqrt(normLib) * MathF.Sqrt(normObs)), 0f, 1f);
        }

        // ── Sort helper ───────────────────────────────────────────────────────

        private static void SortSpan(Span<float> span)
        {
            for (int i = 1; i < span.Length; i++)
            {
                float key = span[i];
                int j = i - 1;
                while (j >= 0 && span[j] > key) { span[j + 1] = span[j]; j--; }
                span[j + 1] = key;
            }
        }
    }
}