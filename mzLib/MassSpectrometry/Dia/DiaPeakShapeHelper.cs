// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Computes peak shape ratio features from the time×fragment intensity matrix.
    /// Phase 16A, Prompt 1: BoundarySignalRatio [27] and ApexToMeanRatio [28] confirmed
    /// active in the classifier. Both were computed previously but not forwarded.
    /// 
    /// These features measure how "peak-like" the chromatographic signal is:
    /// - BoundarySignalRatio: low values indicate sharp peak boundaries (good)
    /// - ApexToMeanRatio: high values indicate a prominent apex relative to background (good)
    /// 
    /// Called from AssembleResultsWithTemporalScoring() when a valid peak group is detected.
    /// </summary>
    internal static class DiaPeakShapeHelper
    {
        /// <summary>
        /// Computes peak shape ratio features: BoundarySignalRatio and ApexToMeanRatio.
        /// Only meaningful when a valid peak group was detected.
        /// 
        /// BoundarySignalRatio = average boundary TIC / apex TIC.
        ///   For true chromatographic peaks, boundaries have much less signal than apex (ratio &lt;&lt; 1).
        ///   For interference or noise, signal is relatively flat (ratio ≈ 1).
        ///   Lower is better for true peptides.
        /// 
        /// ApexToMeanRatio = apex TIC / mean TIC across peak range.
        ///   Higher means a more prominent, sharper peak relative to background.
        ///   Higher is better for true peptides.
        /// </summary>
        /// <param name="matrix">Flat row-major intensity matrix [timePointCount × fragmentCount].</param>
        /// <param name="fragmentCount">Number of fragment columns in the matrix.</param>
        /// <param name="peakLeftIndex">Left boundary scan index of the detected peak group (inclusive).</param>
        /// <param name="peakRightIndex">Right boundary scan index of the detected peak group (inclusive).</param>
        /// <param name="apexIndex">Apex scan index within the peak group.</param>
        /// <param name="result">DiaSearchResult to populate with BoundarySignalRatio and ApexToMeanRatio.</param>
        public static void ComputePeakShapeRatios(
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int peakLeftIndex,
            int peakRightIndex,
            int apexIndex,
            DiaSearchResult result)
        {
            // Compute TIC at the apex scan
            float apexTic = ComputeRowTic(matrix, apexIndex, fragmentCount);

            if (apexTic <= 0f)
                return; // Leave both features as NaN

            // ── BoundarySignalRatio ─────────────────────────────────────────
            float leftTic = ComputeRowTic(matrix, peakLeftIndex, fragmentCount);
            float rightTic = ComputeRowTic(matrix, peakRightIndex, fragmentCount);
            float boundaryTic = (leftTic + rightTic) * 0.5f;

            result.BoundarySignalRatio = boundaryTic / apexTic;

            // ── ApexToMeanRatio ─────────────────────────────────────────────
            // Mean TIC across all scans in [peakLeftIndex, peakRightIndex] inclusive
            int scanCount = peakRightIndex - peakLeftIndex + 1;

            if (scanCount <= 0)
                return; // Degenerate peak group; leave ApexToMeanRatio as NaN

            float sumTic = 0f;
            for (int s = peakLeftIndex; s <= peakRightIndex; s++)
            {
                sumTic += ComputeRowTic(matrix, s, fragmentCount);
            }

            float meanTic = sumTic / scanCount;

            if (meanTic > 0f)
            {
                result.ApexToMeanRatio = apexTic / meanTic;
            }
            // else: leave ApexToMeanRatio as NaN (all-zero peak region)
        }

        /// <summary>
        /// Computes TIC (total ion current = sum of all fragment intensities) for a single
        /// row (scan) in the flat row-major intensity matrix.
        /// </summary>
        /// <param name="matrix">Flat row-major intensity matrix.</param>
        /// <param name="rowIndex">The scan/row index.</param>
        /// <param name="fragmentCount">Number of fragment columns per row.</param>
        /// <returns>Sum of intensities across all fragments for the given row.</returns>
        private static float ComputeRowTic(ReadOnlySpan<float> matrix, int rowIndex, int fragmentCount)
        {
            int offset = rowIndex * fragmentCount;
            ReadOnlySpan<float> row = matrix.Slice(offset, fragmentCount);

            float tic = 0f;
            for (int f = 0; f < row.Length; f++)
            {
                tic += row[f];
            }
            return tic;
        }
    }
}