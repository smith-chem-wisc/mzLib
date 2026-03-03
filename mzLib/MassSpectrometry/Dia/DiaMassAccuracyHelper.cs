// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Phase 13, Prompt 3: Mass accuracy features at the apex scan.
    /// 
    /// For each fragment ion, performs a binary search in the apex scan's m/z array
    /// to find the closest observed peak within the ppm tolerance window.
    /// Computes signed ppm errors and summary statistics across all detected fragments.
    /// 
    /// True peptide signals cluster near the instrument calibration offset (typically near 0 ppm).
    /// Decoy matches scatter randomly across the full tolerance window.
    /// </summary>
    internal static class DiaMassAccuracyHelper
    {
        /// <summary>
        /// Computes mass accuracy features at the apex scan and populates the result.
        /// 
        /// Algorithm:
        ///   1. Get the apex scan's m/z array from DiaScanIndex
        ///   2. For each library fragment m/z, binary search for the closest peak
        ///   3. If within ppm tolerance, record the signed ppm error
        ///   4. Compute mean, median, std, max-abs of the ppm errors
        ///   5. Store observed m/z values for each fragment
        /// 
        /// Requires ≥2 detected fragments at apex for meaningful statistics.
        /// </summary>
        /// <param name="result">DiaSearchResult to populate.</param>
        /// <param name="input">Library precursor input (provides fragment m/z values).</param>
        /// <param name="index">DIA scan index (provides scan-level m/z data).</param>
        /// <param name="apexScanGlobalIndex">
        /// Global scan index in DiaScanIndex for the apex scan.
        /// Must be ≥0 (caller should verify before calling).
        /// </param>
        /// <param name="ppmTolerance">PPM tolerance for matching observed to library m/z.</param>
        public static void ComputeMassAccuracyAtApex(
            DiaSearchResult result,
            in LibraryPrecursorInput input,
            DiaScanIndex index,
            int apexScanGlobalIndex,
            float ppmTolerance)
        {
            if (apexScanGlobalIndex < 0)
                return;

            ReadOnlySpan<float> scanMzs = index.GetScanMzSpan(apexScanGlobalIndex);
            if (scanMzs.Length == 0)
                return;

            int fragmentCount = input.FragmentCount;
            float[] observedMzs = new float[fragmentCount];
            Span<float> ppmErrors = stackalloc float[fragmentCount];
            int nDetected = 0;

            for (int f = 0; f < fragmentCount; f++)
            {
                observedMzs[f] = float.NaN;

                float targetMz = input.FragmentMzs[f];
                float tolDa = targetMz * ppmTolerance * 1e-6f;
                float mzLow = targetMz - tolDa;
                float mzHigh = targetMz + tolDa;

                // Binary search for the insertion point of mzLow
                int lo = 0, hi = scanMzs.Length - 1;
                while (lo <= hi)
                {
                    int mid = lo + (hi - lo) / 2;
                    if (scanMzs[mid] < mzLow)
                        lo = mid + 1;
                    else
                        hi = mid - 1;
                }

                // Search from lo forward for the closest peak within tolerance
                float bestDist = float.MaxValue;
                int bestIdx = -1;
                for (int i = lo; i < scanMzs.Length && scanMzs[i] <= mzHigh; i++)
                {
                    float dist = MathF.Abs(scanMzs[i] - targetMz);
                    if (dist < bestDist)
                    {
                        bestDist = dist;
                        bestIdx = i;
                    }
                }

                if (bestIdx >= 0)
                {
                    float observedMz = scanMzs[bestIdx];
                    observedMzs[f] = observedMz;
                    float ppmError = (observedMz - targetMz) / targetMz * 1e6f;
                    ppmErrors[nDetected++] = ppmError;
                }
            }

            result.ApexObservedMzs = observedMzs;

            if (nDetected < 2)
                return;

            // Compute statistics over the detected ppm errors
            Span<float> detected = ppmErrors.Slice(0, nDetected);

            // Mean
            float sum = 0f;
            for (int i = 0; i < nDetected; i++)
                sum += detected[i];
            float mean = sum / nDetected;
            result.MeanMassErrorPpm = mean;

            // Std (population) and MaxAbs
            float sumSqDev = 0f;
            float maxAbs = 0f;
            for (int i = 0; i < nDetected; i++)
            {
                float dev = detected[i] - mean;
                sumSqDev += dev * dev;
                float absVal = MathF.Abs(detected[i]);
                if (absVal > maxAbs) maxAbs = absVal;
            }
            result.MassErrorStdPpm = MathF.Sqrt(sumSqDev / nDetected);
            result.MaxAbsMassErrorPpm = maxAbs;

            // Median (sort a copy)
            SortSpan(detected);
            result.MedianMassErrorPpm = nDetected % 2 == 1
                ? detected[nDetected / 2]
                : (detected[nDetected / 2 - 1] + detected[nDetected / 2]) / 2f;
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