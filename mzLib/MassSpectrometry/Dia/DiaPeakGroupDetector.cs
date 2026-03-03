// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
// Location: MassSpectrometry/Dia/Scoring/DiaPeakGroupDetector.cs

using System;
using System.Buffers;
using System.Runtime.CompilerServices;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Detected chromatographic peak group for a DIA precursor.
    /// Contains the refined apex, boundaries, and shape features.
    /// </summary>
    public readonly struct PeakGroup
    {
        /// <summary>Index into the time-point array where the apex was found.</summary>
        public readonly int ApexIndex;

        /// <summary>Retention time at the apex (minutes).</summary>
        public readonly float ApexRt;

        /// <summary>Index of the left boundary time point (inclusive).</summary>
        public readonly int LeftIndex;

        /// <summary>Index of the right boundary time point (inclusive).</summary>
        public readonly int RightIndex;

        /// <summary>RT at the left boundary (minutes).</summary>
        public readonly float LeftRt;

        /// <summary>RT at the right boundary (minutes).</summary>
        public readonly float RightRt;

        /// <summary>Peak width in minutes (RightRt - LeftRt).</summary>
        public readonly float PeakWidthMinutes;

        /// <summary>
        /// Symmetry ratio: (apex - left) / (right - left).
        /// 0.5 = perfectly symmetric. &lt;0.5 = left-skewed. &gt;0.5 = right-skewed.
        /// </summary>
        public readonly float SymmetryRatio;

        /// <summary>Number of scans within the peak boundaries (inclusive).</summary>
        public readonly int ScanCount;

        /// <summary>
        /// Total signal (sum of composite intensities) within the peak boundaries.
        /// Higher = more confident detection.
        /// </summary>
        public readonly float TotalSignal;

        /// <summary>
        /// Number of candidate peaks detected before selecting this one.
        /// More candidates = more ambiguous identification.
        /// </summary>
        public readonly int CandidateCount;

        /// <summary>Whether this peak group represents valid detection.</summary>
        public readonly bool IsValid;

        /// <summary>Sentinel for no peak detected.</summary>
        public static readonly PeakGroup None = new PeakGroup();

        public PeakGroup(
            int apexIndex, float apexRt,
            int leftIndex, int rightIndex,
            float leftRt, float rightRt,
            float peakWidthMinutes, float symmetryRatio,
            int scanCount, float totalSignal,
            int candidateCount)
        {
            ApexIndex = apexIndex;
            ApexRt = apexRt;
            LeftIndex = leftIndex;
            RightIndex = rightIndex;
            LeftRt = leftRt;
            RightRt = rightRt;
            PeakWidthMinutes = peakWidthMinutes;
            SymmetryRatio = symmetryRatio;
            ScanCount = scanCount;
            TotalSignal = totalSignal;
            CandidateCount = candidateCount;
            IsValid = true;
        }
    }

    /// <summary>
    /// Detects chromatographic peak groups from DIA precursor XIC data.
    /// 
    /// Algorithm:
    ///   1. Build a composite XIC by summing library-weighted fragment intensities at each time point
    ///   2. Smooth with a moving-average kernel to reduce noise
    ///   3. Find local maxima as candidate peak apices
    ///   4. For each candidate, find boundaries by descending to a fraction of apex height
    ///   5. Score candidates by total signal × scan count, select the best
    /// 
    /// The detected peak boundaries are then used by downstream scoring to restrict
    /// temporal cosine and fragment correlation computations to the actual elution peak,
    /// dramatically reducing interference from flanking co-eluting species.
    /// 
    /// Performance: O(T × F) per precursor where T = time points, F = fragments.
    /// No heap allocations in the detection loop (uses stackalloc for small arrays,
    /// ArrayPool for larger ones).
    /// 
    /// Thread-safe: stateless, operates only on input spans.
    /// </summary>
    public static class DiaPeakGroupDetector
    {
        /// <summary>
        /// Minimum fraction of apex height for peak boundary detection.
        /// A boundary is placed where the smoothed composite drops below this fraction of apex.
        /// 0.05 = 5% of apex height (generous, captures full peak).
        /// </summary>
        private const float BoundaryFraction = 0.05f;

        /// <summary>
        /// Minimum fraction of the global maximum for a local max to be considered a candidate peak.
        /// Prevents noise spikes from being treated as real peaks.
        /// </summary>
        private const float MinCandidateFraction = 0.10f;

        /// <summary>
        /// Minimum number of scans a peak must span to be considered valid.
        /// Peaks narrower than this are likely noise or interference spikes.
        /// </summary>
        private const int MinPeakScans = 3;

        /// <summary>
        /// Half-width of the moving-average smoothing kernel (in time points).
        /// Total kernel width = 2 * SmoothHalfWidth + 1.
        /// 2 = 5-point moving average, suitable for ~14–20 time points per window.
        /// </summary>
        private const int SmoothHalfWidth = 2;

        /// <summary>
        /// Maximum number of candidate peaks to evaluate.
        /// In practice, a ±1 min window rarely has more than 3–4 meaningful peaks.
        /// </summary>
        private const int MaxCandidates = 8;

        /// <summary>
        /// Detects the best peak group from a precursor's time × fragment intensity matrix.
        /// 
        /// The matrix is in row-major order: matrix[t * fragmentCount + f] = intensity
        /// of fragment f at time point t. This is the same matrix built by DiaTemporalScorer.
        /// 
        /// The reference RTs correspond to the time points (rows) of the matrix.
        /// Library intensities are used for weighting the composite XIC.
        /// </summary>
        /// <param name="matrix">Time × fragment intensity matrix, row-major.</param>
        /// <param name="refRts">RT values for each time point (row).</param>
        /// <param name="libraryIntensities">Library intensities per fragment for weighting.</param>
        /// <param name="fragmentCount">Number of fragments (columns).</param>
        /// <param name="timePointCount">Number of time points (rows).</param>
        /// <returns>The best detected peak group, or PeakGroup.None if no valid peak found.</returns>
        public static PeakGroup Detect(
            ReadOnlySpan<float> matrix,
            ReadOnlySpan<float> refRts,
            ReadOnlySpan<float> libraryIntensities,
            int fragmentCount,
            int timePointCount)
        {
            if (timePointCount < MinPeakScans || fragmentCount < 2)
                return PeakGroup.None;

            // ── Step 1: Build composite XIC ─────────────────────────────────
            // Library-weighted sum of fragment intensities at each time point.
            // Using library intensities as weights emphasizes the expected strongest
            // fragments, making the composite more representative of the target.
            bool useStackAlloc = timePointCount <= 128;
            float[] compositeRented = useStackAlloc ? null : ArrayPool<float>.Shared.Rent(timePointCount);
            Span<float> composite = useStackAlloc
                ? stackalloc float[timePointCount]
                : compositeRented.AsSpan(0, timePointCount);

            float[] smoothedRented = useStackAlloc ? null : ArrayPool<float>.Shared.Rent(timePointCount);
            Span<float> smoothed = useStackAlloc
                ? stackalloc float[timePointCount]
                : smoothedRented.AsSpan(0, timePointCount);

            try
            {
                // Compute L2-normalized library weights
                float libNorm = 0f;
                for (int f = 0; f < fragmentCount; f++)
                {
                    float lib = f < libraryIntensities.Length ? libraryIntensities[f] : 0f;
                    libNorm += lib * lib;
                }
                libNorm = libNorm > 0f ? MathF.Sqrt(libNorm) : 1f;

                for (int t = 0; t < timePointCount; t++)
                {
                    float sum = 0f;
                    int rowOffset = t * fragmentCount;
                    for (int f = 0; f < fragmentCount; f++)
                    {
                        float obs = matrix[rowOffset + f];
                        if (obs > 0f)
                        {
                            float libWeight = f < libraryIntensities.Length
                                ? libraryIntensities[f] / libNorm
                                : 1f / fragmentCount;
                            sum += obs * libWeight;
                        }
                    }
                    composite[t] = sum;
                }

                // ── Step 2: Smooth ──────────────────────────────────────────
                MovingAverageSmooth(composite, smoothed, timePointCount, SmoothHalfWidth);

                // ── Step 3: Find candidate peaks ────────────────────────────
                float globalMax = 0f;
                for (int t = 0; t < timePointCount; t++)
                    if (smoothed[t] > globalMax) globalMax = smoothed[t];

                if (globalMax <= 0f)
                    return PeakGroup.None;

                float minCandidateHeight = globalMax * MinCandidateFraction;

                // Stack-allocate candidate storage (MaxCandidates is small)
                Span<int> candidateApices = stackalloc int[MaxCandidates];
                int candidateCount = 0;

                for (int t = 1; t < timePointCount - 1 && candidateCount < MaxCandidates; t++)
                {
                    if (smoothed[t] >= minCandidateHeight &&
                        smoothed[t] >= smoothed[t - 1] &&
                        smoothed[t] >= smoothed[t + 1])
                    {
                        // Ensure it's a true local maximum (not a plateau interior)
                        if (smoothed[t] > smoothed[t - 1] || smoothed[t] > smoothed[t + 1])
                        {
                            candidateApices[candidateCount++] = t;
                        }
                    }
                }

                // Edge case: if no interior local max found, use the global max point
                if (candidateCount == 0)
                {
                    int maxIdx = 0;
                    for (int t = 1; t < timePointCount; t++)
                        if (smoothed[t] > smoothed[maxIdx]) maxIdx = t;

                    if (smoothed[maxIdx] >= minCandidateHeight)
                    {
                        candidateApices[0] = maxIdx;
                        candidateCount = 1;
                    }
                    else
                    {
                        return PeakGroup.None;
                    }
                }

                // ── Step 4: Find boundaries and score each candidate ────────
                PeakGroup bestPeak = PeakGroup.None;
                float bestScore = float.MinValue;

                for (int c = 0; c < candidateCount; c++)
                {
                    int apexIdx = candidateApices[c];
                    float apexHeight = smoothed[apexIdx];
                    float boundaryThreshold = apexHeight * BoundaryFraction;

                    // Find left boundary: descend from apex to left
                    int leftIdx = apexIdx;
                    for (int t = apexIdx - 1; t >= 0; t--)
                    {
                        if (smoothed[t] < boundaryThreshold)
                        {
                            leftIdx = t + 1;
                            break;
                        }
                        // Also stop at a local minimum (valley between peaks)
                        if (t > 0 && smoothed[t] <= smoothed[t - 1] && smoothed[t] < smoothed[t + 1])
                        {
                            leftIdx = t;
                            break;
                        }
                        leftIdx = t;
                    }

                    // Find right boundary: descend from apex to right
                    int rightIdx = apexIdx;
                    for (int t = apexIdx + 1; t < timePointCount; t++)
                    {
                        if (smoothed[t] < boundaryThreshold)
                        {
                            rightIdx = t - 1;
                            break;
                        }
                        // Also stop at a local minimum
                        if (t < timePointCount - 1 && smoothed[t] <= smoothed[t + 1] && smoothed[t] < smoothed[t - 1])
                        {
                            rightIdx = t;
                            break;
                        }
                        rightIdx = t;
                    }

                    int scanCount = rightIdx - leftIdx + 1;
                    if (scanCount < MinPeakScans)
                        continue;

                    // Compute total signal within boundaries (from raw composite, not smoothed)
                    float totalSignal = 0f;
                    for (int t = leftIdx; t <= rightIdx; t++)
                        totalSignal += composite[t];

                    // Peak scoring: prefer high signal with reasonable width
                    float peakScore = totalSignal;

                    if (peakScore > bestScore)
                    {
                        bestScore = peakScore;

                        float apexRt = refRts[apexIdx];
                        float leftRt = refRts[leftIdx];
                        float rightRt = refRts[rightIdx];
                        float width = rightRt - leftRt;
                        float symmetry = width > 0f
                            ? (apexRt - leftRt) / width
                            : 0.5f;

                        bestPeak = new PeakGroup(
                            apexIndex: apexIdx,
                            apexRt: apexRt,
                            leftIndex: leftIdx,
                            rightIndex: rightIdx,
                            leftRt: leftRt,
                            rightRt: rightRt,
                            peakWidthMinutes: width,
                            symmetryRatio: symmetry,
                            scanCount: scanCount,
                            totalSignal: totalSignal,
                            candidateCount: candidateCount);
                    }
                }

                return bestPeak;
            }
            finally
            {
                if (compositeRented != null) ArrayPool<float>.Shared.Return(compositeRented);
                if (smoothedRented != null) ArrayPool<float>.Shared.Return(smoothedRented);
            }
        }

        /// <summary>
        /// Moving-average smoothing with edge-aware kernel.
        /// At boundaries, the kernel is truncated (only valid neighbors contribute).
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void MovingAverageSmooth(
            ReadOnlySpan<float> input, Span<float> output,
            int length, int halfWidth)
        {
            for (int t = 0; t < length; t++)
            {
                int lo = Math.Max(0, t - halfWidth);
                int hi = Math.Min(length - 1, t + halfWidth);
                float sum = 0f;
                int count = hi - lo + 1;
                for (int k = lo; k <= hi; k++)
                    sum += input[k];
                output[t] = sum / count;
            }
        }
    }
}