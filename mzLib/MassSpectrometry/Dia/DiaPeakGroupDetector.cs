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
    /// 
    /// Prompt 2 additions:
    ///   - ApexCosine: cosine similarity between library and observed at the apex
    ///   - CoElutionStd: std dev of per-fragment apex positions within peak boundaries
    ///   - Snr: composite apex / flanking baseline ratio
    ///   - SelectionScore: ApexCosine × CoElutionFactor × SnrFactor
    ///   - TotalCandidateCount: how many candidates existed before selection
    ///   - SecondBestScore: SelectionScore of the runner-up (0 if only 1 candidate)
    /// 
    /// All fields from the original PeakGroup struct are preserved for backward compatibility
    /// with DiaFeatureExtractor and DiaSearchResult callers.
    /// </summary>
    public readonly struct PeakGroup
    {
        // ════════════════════════════════════════════════════════════════
        //  Original PeakGroup fields — preserved with identical semantics
        // ════════════════════════════════════════════════════════════════

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
        /// RT of every candidate peak apex found during detection, before selection.
        /// Used by grid calibration to take median over all candidates across a bin
        /// rather than committing to one apex per precursor.
        /// </summary>
        public readonly float[] CandidateApexRts;

        /// <summary>
        /// Number of candidate peaks detected before selecting this one.
        /// More candidates = more ambiguous identification.
        /// </summary>
        public readonly int CandidateCount;

        /// <summary>Whether this peak group represents valid detection.</summary>
        public readonly bool IsValid;

        // ════════════════════════════════════════════════════════════════
        //  NEW fields (Prompt 2)
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Cosine similarity between library fragment intensities and observed fragment
        /// intensities at the apex time point.
        /// Range [0, 1]. Higher = better spectral match at apex.
        /// </summary>
        public readonly float ApexCosine;

        /// <summary>
        /// Standard deviation (in scan units) of per-fragment apex positions within
        /// the peak boundaries. Low value = all fragments co-elute tightly (true peptide peak).
        /// High value = fragments peak at different times (interference).
        /// </summary>
        public readonly float CoElutionStd;

        /// <summary>
        /// Signal-to-noise ratio at the apex: composite apex / flanking baseline.
        /// Higher = cleaner, more prominent peak.
        /// </summary>
        public readonly float Snr;

        /// <summary>
        /// Composite selection score: ApexCosine × CoElutionFactor × SnrFactor.
        /// This is the score used to pick the best candidate among all detected peaks.
        /// Range [0, 1].
        /// </summary>
        public readonly float SelectionScore;

        /// <summary>
        /// Total number of candidate peaks found before selection.
        /// (Same as CandidateCount but kept as a separate field per spec.)
        /// </summary>
        public readonly int TotalCandidateCount;

        /// <summary>
        /// SelectionScore of the runner-up candidate (0 if only one candidate existed).
        /// Large gap (SelectionScore - SecondBestScore) indicates confident selection.
        /// </summary>
        public readonly float SecondBestScore;

        // ════════════════════════════════════════════════════════════════
        //  Sentinel
        // ════════════════════════════════════════════════════════════════

        /// <summary>Sentinel for no peak detected.</summary>
        public static readonly PeakGroup None = new PeakGroup();

        // ════════════════════════════════════════════════════════════════
        //  Constructors
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Full constructor used by DiaPeakGroupDetector.SelectBest().
        /// </summary>
        public PeakGroup(
            int apexIndex, float apexRt,
            int leftIndex, int rightIndex,
            float leftRt, float rightRt,
            float peakWidthMinutes, float symmetryRatio,
            int scanCount, float totalSignal,
            int candidateCount,
            float[] candidateApexRts,
            float apexCosine,
            float coElutionStd,
            float snr,
            float selectionScore,
            int totalCandidateCount,
            float secondBestScore)
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
            CandidateApexRts = candidateApexRts;
            ApexCosine = apexCosine;
            CoElutionStd = coElutionStd;
            Snr = snr;
            SelectionScore = selectionScore;
            TotalCandidateCount = totalCandidateCount;
            SecondBestScore = secondBestScore;
            IsValid = true;
        }
    }

    // Backward-compatibility alias — callers that still reference PeakGroup continue to compile.
    // Remove after all callers have been updated to use PeakGroup directly.

    /// <summary>
    /// Detects chromatographic peak groups from DIA precursor XIC data.
    ///
    /// Prompt 2 rewrite: SelectBest() replaces Detect() as the primary entry point.
    ///
    /// New algorithm:
    ///   1. Build library-weighted composite XIC and smooth (same as before)
    ///   2. Find candidate local maxima above MinCandidateFraction threshold (same)
    ///   3. For each candidate, score on THREE axes:
    ///      - ApexCosine: spectral match at the apex time point
    ///      - CoElutionFactor: exp(-coElutionStd / CoElutionSigmaScans) — penalizes
    ///        candidates where fragment XICs peak at different time points
    ///      - SnrFactor: log2(1 + snr) / log2(1 + SnrNormalizationTarget) — rewards
    ///        candidates with high apex signal relative to baseline
    ///   4. SelectionScore = ApexCosine × CoElutionFactor × SnrFactor
    ///   5. Pick the candidate with highest SelectionScore (RT tiebreak optional)
    ///
    /// Why this fixes the "4.7 min early" systematic error:
    ///   The old scorer (cosine × log(signal)) selected interference peaks that had
    ///   both high intensity AND high cosine at the wrong RT. The new scorer adds
    ///   CoElutionFactor which penalizes any candidate where the individual fragment
    ///   XICs do not agree on a single elution apex. Interference peaks tend to pull
    ///   different fragments to different time points, producing a high CoElutionStd
    ///   → low CoElutionFactor → lower SelectionScore despite high cosine.
    ///
    /// Performance: O(T × F × C) per precursor where C ≤ MaxCandidates = 20.
    /// Thread-safe: stateless, operates only on input spans.
    /// </summary>
    public static class DiaPeakGroupDetector
    {
        // ════════════════════════════════════════════════════════════════
        //  Configuration constants
        // ════════════════════════════════════════════════════════════════

        private const int MaxCandidates = 20;   // raised from 5; wide windows need more candidates
        private const float MinCandidateFraction = 0.10f;
        private const int SmoothHalfWidth = 2;
        private const float BoundaryFraction = 0.05f;
        private const float CoElutionSigmaScans = 1.5f;
        private const float SnrNormalizationTarget = 10.0f;
        private const int MinPeakScans = 3;

        // ════════════════════════════════════════════════════════════════
        //  Primary entry point — Prompt 2 multi-hypothesis scorer
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Selects the best peak group from a precursor's time × fragment intensity matrix
        /// using a three-axis scoring model: spectral match × co-elution consistency × SNR.
        ///
        /// The matrix is in row-major order: matrix[t * fragmentCount + f] = intensity
        /// of fragment f at time point t.
        /// </summary>
        /// <param name="matrix">Time × fragment intensity matrix, row-major.</param>
        /// <param name="refRts">RT values for each time point (row), in minutes.</param>
        /// <param name="libraryIntensities">Library fragment intensities for weighting/scoring.</param>
        /// <param name="fragmentCount">Number of fragments (columns).</param>
        /// <param name="timePointCount">Number of time points (rows).</param>
        /// <param name="predictedRt">
        /// Predicted retention time from the RT calibration model (minutes).
        /// When provided, the RT proximity Gaussian is the PRIMARY discriminating factor
        /// between candidates — not a tiebreaker. The Gaussian weight crushes candidates
        /// far from this RT (e.g. a peak 4.7 min away gets a ~62,000× penalty relative
        /// to the true peak 0.1 min away at ±5 min window).
        /// Passing null disables RT-based scoring entirely, which is only appropriate
        /// for libraries with no RT annotation.
        /// </param>
        /// <returns>
        /// The best detected PeakGroup, or PeakGroup.None if no valid
        /// peak was found (too few time points or no signal).
        /// </returns>
        public static PeakGroup SelectBest(
            ReadOnlySpan<float> matrix,
            ReadOnlySpan<float> refRts,
            ReadOnlySpan<float> libraryIntensities,
            int fragmentCount,
            int timePointCount,
            float? predictedRt = null,
            float rtWindowHalfWidth = 0f)
        {
            if (timePointCount < MinPeakScans || fragmentCount < 2)
                return PeakGroup.None;

            // ── Step 1: Build library-weighted composite XIC ─────────────
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
                // Library-weighted sum: weight by sqrt(library intensity) for composite XIC.
                // This emphasises the library's dominant fragments without over-weighting them.
                for (int t = 0; t < timePointCount; t++)
                {
                    float sum = 0f;
                    int rowOffset = t * fragmentCount;
                    for (int f = 0; f < fragmentCount; f++)
                    {
                        float obs = matrix[rowOffset + f];
                        if (obs > 0f)
                        {
                            float lib = f < libraryIntensities.Length ? libraryIntensities[f] : 0f;
                            sum += obs * MathF.Sqrt(lib);
                        }
                    }
                    composite[t] = sum;
                }

                // ── Step 1b: Smooth ───────────────────────────────────────
                MovingAverageSmooth(composite, smoothed, timePointCount, SmoothHalfWidth);

                // ── Step 2: Find candidate peaks ──────────────────────────
                float globalMax = 0f;
                for (int t = 0; t < timePointCount; t++)
                    if (smoothed[t] > globalMax) globalMax = smoothed[t];

                if (globalMax <= 0f)
                    return PeakGroup.None;

                float minCandidateHeight = globalMax * MinCandidateFraction;

                // Fixed-size stack buffer for candidate indices
                Span<int> candidateApices = stackalloc int[MaxCandidates];
                int candidateCount = 0;

                for (int t = 1; t < timePointCount - 1 && candidateCount < MaxCandidates; t++)
                {
                    if (smoothed[t] >= minCandidateHeight &&
                        smoothed[t] >= smoothed[t - 1] &&
                        smoothed[t] >= smoothed[t + 1])
                    {
                        if (smoothed[t] > smoothed[t - 1] || smoothed[t] > smoothed[t + 1])
                            candidateApices[candidateCount++] = t;
                    }
                }

                // Edge case: no interior local max — use global max
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

                // Collect all candidate apex RTs (for CandidateApexRts field, used by calibration)
                float[] candidateApexRts = new float[candidateCount];
                for (int c = 0; c < candidateCount; c++)
                    candidateApexRts[c] = refRts[candidateApices[c]];

                // ── Step 3: Pre-compute SNR baseline once ─────────────────
                // Baseline = mean composite signal in outer 20% of the FULL extraction window
                // (not per-peak — using global flanking to avoid peak contaminating baseline)
                int outerBand = Math.Max(1, timePointCount / 5);
                float baselineSum = 0f;
                int baselineCount = 0;
                for (int t = 0; t < outerBand; t++)
                {
                    baselineSum += composite[t];
                    baselineCount++;
                }
                for (int t = timePointCount - outerBand; t < timePointCount; t++)
                {
                    baselineSum += composite[t];
                    baselineCount++;
                }
                float baseline = baselineCount > 0 ? baselineSum / baselineCount : 0f;
                if (baseline < 1e-6f) baseline = 1e-6f;

                // ── Step 4: Score each candidate ──────────────────────────
                // Per-candidate score storage (all value types — safe for stackalloc)
                Span<int> candApex = stackalloc int[MaxCandidates];
                Span<int> candLeft = stackalloc int[MaxCandidates];
                Span<int> candRight = stackalloc int[MaxCandidates];
                Span<float> candScore = stackalloc float[MaxCandidates];
                Span<float> candCosine = stackalloc float[MaxCandidates];
                Span<float> candCoElut = stackalloc float[MaxCandidates];
                Span<float> candSnr = stackalloc float[MaxCandidates];
                int validCandCount = 0;

                for (int c = 0; c < candidateCount; c++)
                {
                    int apexIdx = candidateApices[c];
                    float apexHeight = smoothed[apexIdx];
                    float boundaryThr = apexHeight * BoundaryFraction;

                    // Find left boundary
                    int leftIdx = apexIdx;
                    for (int t = apexIdx - 1; t >= 0; t--)
                    {
                        if (smoothed[t] < boundaryThr)
                        {
                            leftIdx = t + 1;
                            break;
                        }
                        if (t > 0 && smoothed[t] <= smoothed[t - 1] && smoothed[t] < smoothed[t + 1])
                        {
                            leftIdx = t;
                            break;
                        }
                        leftIdx = t;
                    }

                    // Find right boundary
                    int rightIdx = apexIdx;
                    for (int t = apexIdx + 1; t < timePointCount; t++)
                    {
                        if (smoothed[t] < boundaryThr)
                        {
                            rightIdx = t - 1;
                            break;
                        }
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

                    // ── Score 1: ApexCosine ──────────────────────────────
                    int apexRowOffset = apexIdx * fragmentCount;
                    float dot = 0f, normLib = 0f, normObs = 0f;
                    int cosineContributors = 0;
                    for (int f = 0; f < fragmentCount; f++)
                    {
                        float obs = matrix[apexRowOffset + f];
                        if (obs <= 0f) continue;
                        float lib = f < libraryIntensities.Length ? libraryIntensities[f] : 0f;
                        if (lib <= 0f) continue;
                        dot += lib * obs;
                        normLib += lib * lib;
                        normObs += obs * obs;
                        cosineContributors++;
                    }
                    float apexCosine = (cosineContributors >= 2 && normLib > 0f && normObs > 0f)
                        ? Math.Clamp(dot / (MathF.Sqrt(normLib) * MathF.Sqrt(normObs)), 0f, 1f)
                        : 0f;

                    // ── Score 2: CoElution consistency ───────────────────
                    // For each fragment, find the time index of its maximum intensity
                    // within [leftIdx, rightIdx].
                    // CoElutionStd = std dev of those per-fragment apex indices.
                    int coElutFragCount = 0;
                    float coElutSum = 0f;
                    float coElutSumSq = 0f;

                    for (int f = 0; f < fragmentCount; f++)
                    {
                        float maxInt = 0f;
                        int maxT = leftIdx;
                        for (int t = leftIdx; t <= rightIdx; t++)
                        {
                            float v = matrix[t * fragmentCount + f];
                            if (v > maxInt) { maxInt = v; maxT = t; }
                        }
                        if (maxInt > 0f)
                        {
                            float fi = (float)maxT;
                            coElutSum += fi;
                            coElutSumSq += fi * fi;
                            coElutFragCount++;
                        }
                    }

                    float coElutionStd = 0f;
                    if (coElutFragCount >= 2)
                    {
                        float mean = coElutSum / coElutFragCount;
                        float variance = coElutSumSq / coElutFragCount - mean * mean;
                        coElutionStd = variance > 0f ? MathF.Sqrt(variance) : 0f;
                    }
                    float coElutionFactor = MathF.Exp(-coElutionStd / CoElutionSigmaScans);

                    // ── Score 3: SNR ─────────────────────────────────────
                    float snr = composite[apexIdx] / baseline;
                    float snrFactor = MathF.Log2(1f + snr) / MathF.Log2(1f + SnrNormalizationTarget);
                    snrFactor = Math.Clamp(snrFactor, 0f, 1f);

                    // ── Combined SelectionScore ──────────────────────────
                    // Base score: spectral match × co-elution coherence × signal prominence
                    float baseScore = apexCosine * coElutionFactor * snrFactor;

                    // ── RT proximity factor ───────────────────────────────
                    // Gaussian weight centred on predictedRt.
                    // This is the primary discriminator when interference peaks are present:
                    // a peak 4.7 min from the library RT is exponentially crushed relative
                    // to the true peak ~0.1 min from the library RT.
                    //
                    // When predictedRt is null or rtWindowHalfWidth == 0,
                    // rtProximityFactor = 1.0 and RT-based scoring is disabled.
                    float rtProximityFactor = 1.0f;
                    if (predictedRt.HasValue && rtWindowHalfWidth > 0f)
                    {
                        // rtSigma = max(halfWidth * 0.20, 0.5 min)
                        // This makes the RT proximity factor overwhelmingly decisive:
                        //   ±5 min window (sigma=1.0): wrong peak 4.7 min away → factor = 9e-6
                        //   ±1 min window (sigma=0.5): true peak 0.4 min away  → factor = 0.72
                        // The floor of 0.5 min prevents over-penalizing at narrow windows
                        // where calibration may place the iRT 0.3-0.5 min from the true apex.
                        float rtSigma = Math.Max(rtWindowHalfWidth * 0.20f, 0.5f);
                        float deltaRt = MathF.Abs(refRts[apexIdx] - predictedRt.Value);
                        float zScore = deltaRt / rtSigma;
                        rtProximityFactor = MathF.Exp(-0.5f * zScore * zScore);
                    }

                    float selectionScore = baseScore * rtProximityFactor;

                    // Record valid candidate
                    candApex[validCandCount] = apexIdx;
                    candLeft[validCandCount] = leftIdx;
                    candRight[validCandCount] = rightIdx;
                    candScore[validCandCount] = selectionScore;
                    candCosine[validCandCount] = apexCosine;
                    candCoElut[validCandCount] = coElutionStd;
                    candSnr[validCandCount] = snr;
                    validCandCount++;
                }

                if (validCandCount == 0)
                    return PeakGroup.None;

                // ── Step 4: Select winner ─────────────────────────────────
                int bestIdx = 0;
                float bestScore = candScore[0];
                int secondIdx = -1;
                float secondScore = 0f;

                for (int c = 1; c < validCandCount; c++)
                {
                    if (candScore[c] > bestScore)
                    {
                        secondIdx = bestIdx;
                        secondScore = bestScore;
                        bestIdx = c;
                        bestScore = candScore[c];
                    }
                    else if (secondIdx < 0 || candScore[c] > secondScore)
                    {
                        secondIdx = c;
                        secondScore = candScore[c];
                    }
                }

                // ── Step 4b: No separate tiebreak needed ──────────────────
                // RT proximity is now baked into SelectionScore as a continuous
                // Gaussian factor. A separate tiebreak on top would be redundant
                // and could re-introduce the bias it was meant to avoid.

                // ── Step 5: Build the winning PeakGroup ──────────
                int winApex = candApex[bestIdx];
                int winLeft = candLeft[bestIdx];
                int winRight = candRight[bestIdx];
                int winScans = winRight - winLeft + 1;

                float totalSignal = 0f;
                for (int t = winLeft; t <= winRight; t++)
                    totalSignal += composite[t];

                float apexRt = refRts[winApex];
                float leftRt = refRts[winLeft];
                float rightRt = refRts[winRight];
                float width = rightRt - leftRt;
                float symmetry = width > 0f ? (apexRt - leftRt) / width : 0.5f;

                return new PeakGroup(
                    apexIndex: winApex,
                    apexRt: apexRt,
                    leftIndex: winLeft,
                    rightIndex: winRight,
                    leftRt: leftRt,
                    rightRt: rightRt,
                    peakWidthMinutes: width,
                    symmetryRatio: symmetry,
                    scanCount: winScans,
                    totalSignal: totalSignal,
                    candidateCount: candidateCount,
                    candidateApexRts: candidateApexRts,
                    apexCosine: candCosine[bestIdx],
                    coElutionStd: candCoElut[bestIdx],
                    snr: candSnr[bestIdx],
                    selectionScore: bestScore,
                    totalCandidateCount: candidateCount,
                    secondBestScore: secondIdx >= 0 ? secondScore : 0f);
            }
            finally
            {
                if (compositeRented != null) ArrayPool<float>.Shared.Return(compositeRented);
                if (smoothedRented != null) ArrayPool<float>.Shared.Return(smoothedRented);
            }
        }

        // ════════════════════════════════════════════════════════════════
        //  Legacy entry point — preserved for backward compatibility
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Legacy peak detection entry point from Prompt 1.
        /// Now delegates to SelectBest() with no predictedRt (RT scoring disabled).
        /// Preserved so that any callers not yet updated to SelectBest() continue to compile.
        /// </summary>

        public static PeakGroup Detect(
            ReadOnlySpan<float> matrix,
            ReadOnlySpan<float> refRts,
            ReadOnlySpan<float> libraryIntensities,
            int fragmentCount,
            int timePointCount,
            float minCandidateFraction = MinCandidateFraction)
        {
            // minCandidateFraction is ignored here — SelectBest always uses the class constant.
            // The parameter is kept only to avoid breaking old call sites.
            return SelectBest(matrix, refRts, libraryIntensities, fragmentCount, timePointCount,
                predictedRt: null);
        }

        // ════════════════════════════════════════════════════════════════
        //  Shared utilities
        // ════════════════════════════════════════════════════════════════

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