// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Buffers;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Feature vector for a single DIA precursor identification.
    ///
    /// Feature layout history:
    ///   Phase 13 Prompt 8:  26 features [0-25]
    ///   Phase 16A Prompt 1: 29 features — activated BestFragWeightedCosine [26],
    ///                        BoundarySignalRatio [27], ApexToMeanRatio [28]
    ///   Phase 16B Prompt 6: 33 features — MS1 features [29-32]
    ///   Phase 19:           35 features — ChimericScore [33],
    ///                        LibraryCoverageFraction [34]
    ///                        (RtDeviationNormalized dropped: 100% NaN due to
    ///                         PeakWidth=0 when no peak group detected)
    /// </summary>
    public struct DiaFeatureVector
    {
        // ── Primary scores (3) ───────────────────────────────────────
        public float ApexScore;             // [0]
        public float TemporalScore;         // [1]
        public float SpectralAngle;         // [2]

        // ── Peak-region correlations (2) ─────────────────────────────
        public float PeakMeanFragCorr;      // [3]
        public float PeakMinFragCorr;       // [4]

        // ── Peak shape (2) ───────────────────────────────────────────
        public float PeakWidth;             // [5]
        public float CandidateCount;        // [6]

        // ── Signal features (2) ──────────────────────────────────────
        public float LogTotalIntensity;     // [7]
        public float IntensityCV;           // [8]

        // ── Fragment evidence (1) ────────────────────────────────────
        public float FragmentDetectionRate; // [9]

        // ── Retention time (2) ───────────────────────────────────────
        public float RtDeviationMinutes;    // [10]
        public float RtDeviationSquared;    // [11]

        // ── Temporal evidence (1) ────────────────────────────────────
        public int TimePointsUsed;          // [12]

        // ── Mass accuracy (3) ────────────────────────────────────────
        public float MeanMassErrorPpm;      // [13]  (absolute value of mean signed error)
        public float MassErrorStdPpm;       // [14]
        public float MaxAbsMassErrorPpm;    // [15]

        // ── Best-fragment reference curve (4) ────────────────────────
        public float BestFragCorrelationSum;// [16]
        public float MedianFragRefCorr;     // [17]
        public float MinFragRefCorr;        // [18]
        public float StdFragRefCorr;        // [19]

        // ── Signal ratio deviation (3) ───────────────────────────────
        public float MeanSignalRatioDeviation;  // [20]
        public float MaxSignalRatioDeviation;   // [21]
        public float StdSignalRatioDeviation;   // [22]

        // ── Smoothed correlations (2) ────────────────────────────────
        public float SmoothedMeanFragCorr;  // [23]
        public float SmoothedMinFragCorr;   // [24]

        // ── Signal-to-noise (1) ──────────────────────────────────────
        public float Log2SignalToNoise;     // [25]

        // ── Migrated features (3) — Phase 16A ────────────────────────
        public float BestFragWeightedCosine;   // [26]
        public float BoundarySignalRatio;      // [27]
        public float ApexToMeanRatio;          // [28]

        // ── MS1 features (4) — Phase 16B ─────────────────────────────
        public float PrecursorXicApexIntensity; // [29]
        public float IsotopePatternScore;        // [30]
        public float Ms1Ms2Correlation;          // [31]
        public float PrecursorElutionScore;      // [32]

        // ── Interference / chimeric (1) — Phase 19 ───────────────────
        public float ChimericScore;              // [33]

        // ── Coverage (1) — Phase 19 ──────────────────────────────────
        public float LibraryCoverageFraction;    // [34]

        // ── Metadata (not classifier features) ───────────────────────
        public bool IsDecoy;
        public int PrecursorIndex;
        public int ChargeState;
        public float PrecursorMz;
        public int FragmentsDetected;
        public int FragmentsQueried;

        /// <summary>
        /// Number of features used by the classifier.
        /// Phase 13 Prompt 8:  26 features [0-25]
        /// Phase 16A Prompt 1: 29 features (+ migrated [26-28])
        /// Phase 16B Prompt 6: 33 features (+ MS1 [29-32])
        /// Phase 19:           35 features (+ ChimericScore[33] + LibraryCoverageFraction[34])
        /// </summary>
        public const int ClassifierFeatureCount = 35;

        /// <summary>
        /// Index of ApexScore in the feature vector.
        /// Used by DiaLinearDiscriminant for interaction features.
        /// </summary>
        public const int InteractionFeatureIndexA = 0;  // ApexScore

        /// <summary>
        /// Index of PeakMeanFragCorr in the feature vector.
        /// Used by DiaLinearDiscriminant for interaction features.
        /// </summary>
        public const int InteractionFeatureIndexB = 3;  // PeakMeanFragCorr

        /// <summary>
        /// Writes classifier features into a float span.
        /// Order must be consistent with weight vectors.
        ///
        ///   [0-2]   Primary scores: ApexScore, TemporalScore, SpectralAngle
        ///   [3-4]   Peak correlations: PeakMeanFragCorr, PeakMinFragCorr
        ///   [5-6]   Peak shape: PeakWidth, CandidateCount
        ///   [7-8]   Signal: LogTotalIntensity, IntensityCV
        ///   [9]     Fragment evidence: FragDetRate
        ///   [10-11] RT: RtDeviationMinutes, RtDeviationSquared
        ///   [12]    Temporal: TimePointsUsed
        ///   [13-15] Mass accuracy: MeanMassErrorPpm(abs), MassErrorStdPpm, MaxAbsMassErrorPpm
        ///   [16-19] Best-fragment: BestFragCorrSum, MedianRefCorr, MinRefCorr, StdRefCorr
        ///   [20-22] Signal ratio: MeanSigRatioDev, MaxSigRatioDev, StdSigRatioDev
        ///   [23-24] Smoothed correlations: SmoothedMeanFragCorr, SmoothedMinFragCorr
        ///   [25]    S/N: Log2SignalToNoise
        ///   [26-28] Migrated: BestFragWeightedCosine, BoundarySignalRatio, ApexToMeanRatio
        ///   [29-32] MS1: PrecursorXicApexIntensity, IsotopePatternScore,
        ///                Ms1Ms2Correlation, PrecursorElutionScore
        ///   [33]    Interference: ChimericScore
        ///   [34]    Coverage: LibraryCoverageFraction
        /// </summary>
        public readonly void WriteTo(Span<float> features)
        {
            if (features.Length < ClassifierFeatureCount)
                throw new ArgumentException(
                    $"Span must have at least {ClassifierFeatureCount} elements");

            features[0] = ApexScore;
            features[1] = TemporalScore;
            features[2] = SpectralAngle;
            features[3] = PeakMeanFragCorr;
            features[4] = PeakMinFragCorr;
            features[5] = PeakWidth;
            features[6] = CandidateCount;
            features[7] = LogTotalIntensity;
            features[8] = IntensityCV;
            features[9] = FragmentDetectionRate;
            features[10] = RtDeviationMinutes;
            features[11] = RtDeviationSquared;
            features[12] = (float)TimePointsUsed;
            features[13] = MeanMassErrorPpm;
            features[14] = MassErrorStdPpm;
            features[15] = MaxAbsMassErrorPpm;
            features[16] = BestFragCorrelationSum;
            features[17] = MedianFragRefCorr;
            features[18] = MinFragRefCorr;
            features[19] = StdFragRefCorr;
            features[20] = MeanSignalRatioDeviation;
            features[21] = MaxSignalRatioDeviation;
            features[22] = StdSignalRatioDeviation;
            features[23] = SmoothedMeanFragCorr;
            features[24] = SmoothedMinFragCorr;
            features[25] = Log2SignalToNoise;

            // Migrated features [26-28]
            features[26] = BestFragWeightedCosine;
            features[27] = BoundarySignalRatio;
            features[28] = ApexToMeanRatio;

            // MS1 features [29-32]
            features[29] = PrecursorXicApexIntensity;
            features[30] = IsotopePatternScore;
            features[31] = Ms1Ms2Correlation;
            features[32] = PrecursorElutionScore;

            // Interference / chimeric [33]
            features[33] = ChimericScore;

            // Coverage [34]
            features[34] = LibraryCoverageFraction;
        }

        /// <summary>
        /// Feature names for TSV/diagnostic output.
        /// Order matches WriteTo().
        /// </summary>
        public static readonly string[] FeatureNames = new[]
        {
            // Primary scores [0-2]
            "ApexScore", "TemporalScore", "SpectralAngle",
            // Peak correlations [3-4]
            "PeakMeanFragCorr", "PeakMinFragCorr",
            // Peak shape [5-6]
            "PeakWidth", "CandidateCount",
            // Signal [7-8]
            "LogTotalIntensity", "IntensityCV",
            // Fragment evidence [9]
            "FragDetRate",
            // RT [10-11]
            "RtDeviationMinutes", "RtDeviationSquared",
            // Temporal [12]
            "TimePointsUsed",
            // Mass accuracy [13-15]
            "MeanMassErrorPpm", "MassErrorStdPpm", "MaxAbsMassErrorPpm",
            // Best-fragment [16-19]
            "BestFragCorrSum", "MedianFragRefCorr", "MinFragRefCorr", "StdFragRefCorr",
            // Signal ratio [20-22]
            "MeanSigRatioDev", "MaxSigRatioDev", "StdSigRatioDev",
            // Smoothed correlations [23-24]
            "SmoothedMeanFragCorr", "SmoothedMinFragCorr",
            // S/N [25]
            "Log2SNR",
            // Migrated features [26-28]
            "BestFragWeightedCosine", "BoundarySignalRatio", "ApexToMeanRatio",
            // MS1 features [29-32]
            "PrecursorXicApexIntensity", "IsotopePatternScore",
            "Ms1Ms2Correlation", "PrecursorElutionScore",
            // Interference / chimeric [33]
            "ChimericScore",
            // Coverage [34]
            "LibraryCoverageFraction",
        };
    }

    /// <summary>
    /// Computes feature vectors from DiaSearchResult objects.
    /// Thread-safe: uses only local state + ArrayPool.
    ///
    /// Phase 16A, Prompt 1: 29-feature vector with all migrated features active.
    /// Phase 16B, Prompt 6: 33-feature vector. When a DiaScanIndex with MS1 data is
    /// provided, ComputeMs1Features() is called to populate features [29-32] on the
    /// result, then those values are read back into the vector. When index is null or
    /// has no MS1 scans, features [29-32] remain NaN (handled by mean imputation in
    /// the classifier).
    /// Phase 19: 35-feature vector. ChimericScore[33] and LibraryCoverageFraction[34]
    /// are computed during assembly by DiaLibraryQueryGenerator and stored on result.
    /// </summary>
    public static class DiaFeatureExtractor
    {
        /// <summary>
        /// Computes a 35-feature vector from a scored DiaSearchResult.
        ///
        /// MS1 features [29-32] require a DiaScanIndex built with MS1 scans.
        /// Pass null for <paramref name="index"/> to skip MS1 features (they remain NaN).
        /// </summary>
        /// <param name="result">Scored result from the extraction pipeline.</param>
        /// <param name="precursorIndex">Zero-based index of this precursor in the results list.</param>
        /// <param name="index">
        /// Optional DiaScanIndex with MS1 scan data. When provided and Ms1ScanCount > 0,
        /// MS1 features are computed and stored on <paramref name="result"/> before being
        /// read into the vector.
        /// </param>
        /// <param name="bestFragXic">
        /// Optional: intensities of the best MS2 fragment XIC, used for Ms1Ms2Correlation [31].
        /// Length must match <paramref name="bestFragXicRts"/>. Pass default to omit.
        /// </param>
        /// <param name="bestFragXicRts">
        /// Optional: retention times corresponding to <paramref name="bestFragXic"/>.
        /// Pass default to omit.
        /// </param>
        public static DiaFeatureVector ComputeFeatures(
            DiaSearchResult result,
            int precursorIndex,
            DiaScanIndex index = null,
            ReadOnlySpan<float> bestFragXic = default,
            ReadOnlySpan<float> bestFragXicRts = default)
        {
            var fv = new DiaFeatureVector();

            // ── Primary scores [0-2] ────────────────────────────────────────
            fv.ApexScore = SafeScore(result.ApexScore);
            fv.TemporalScore = SafeScore(result.TemporalScore);
            // Coalesce: temporal path sets SpectralAngle; simple assembler sets SpectralAngleScore.
            float _sa = result.SpectralAngle;
            fv.SpectralAngle = SafeScore(float.IsNaN(_sa) ? result.SpectralAngleScore : _sa);

            // ── Peak-region correlations [3-4] ──────────────────────────────
            if (result.DetectedPeakGroup.HasValue && result.DetectedPeakGroup.Value.IsValid)
            {
                fv.PeakMeanFragCorr = SafeScore(result.PeakMeanFragCorr);
                fv.PeakMinFragCorr = float.IsNaN(result.PeakMinFragCorr)
                    ? -1f : result.PeakMinFragCorr;
                fv.PeakWidth = result.DetectedPeakGroup.Value.PeakWidthMinutes;
            }
            else
            {
                // Fallback: use full-window correlations when no peak group detected
                fv.PeakMeanFragCorr = SafeScore(result.MeanFragCorr);
                fv.PeakMinFragCorr = float.IsNaN(result.MinFragCorr)
                    ? -1f : result.MinFragCorr;
                fv.PeakWidth = 0f;
            }

            // ── Peak shape [5-6] ────────────────────────────────────────────
            // PeakWidth already set above
            fv.CandidateCount = SafeScore(result.CandidateCount);

            // ── Fragment evidence [9] ───────────────────────────────────────
            fv.FragmentDetectionRate = result.FragmentDetectionRate;
            fv.FragmentsDetected = result.FragmentsDetected;
            fv.FragmentsQueried = result.FragmentsQueried;

            // ── Temporal evidence [12] ──────────────────────────────────────
            fv.TimePointsUsed = result.TimePointsUsed;

            // ── Intensity features [7-8] ────────────────────────────────────
            float totalIntensity = 0f;
            int nDetected = 0;

            for (int f = 0; f < result.FragmentsQueried; f++)
            {
                totalIntensity += result.ExtractedIntensities[f];
                if (result.XicPointCounts[f] > 0)
                    nDetected++;
            }

            fv.LogTotalIntensity = totalIntensity > 0 ? MathF.Log10(totalIntensity) : 0f;

            if (nDetected >= 2)
            {
                float mean = 0f;
                for (int f = 0; f < result.FragmentsQueried; f++)
                    if (result.XicPointCounts[f] > 0)
                        mean += result.ExtractedIntensities[f];
                mean /= nDetected;

                float variance = 0f;
                for (int f = 0; f < result.FragmentsQueried; f++)
                {
                    if (result.XicPointCounts[f] > 0)
                    {
                        float diff = result.ExtractedIntensities[f] - mean;
                        variance += diff * diff;
                    }
                }
                variance /= nDetected;
                fv.IntensityCV = mean > 0 ? MathF.Sqrt(variance) / mean : 0f;
            }
            else
            {
                fv.IntensityCV = 1f;
            }

            // ── Retention time features [10-11] ─────────────────────────────
            // IMPORTANT: Use the pre-calibrated RtDeviationMinutes set by
            // DiaCalibrationPipeline.RecalibrateRtDeviations(), which converts
            // iRT units → observed minutes via the calibration model.
            // Do NOT recompute from LibraryRetentionTime directly — that value may
            // be in iRT units (range -6 to 121), not observed minutes. Direct
            // subtraction produces deltaRt ~25 min for all results, which clamps
            // to MaxRtDeviationMinutes for every target and decoy, destroying
            // discriminative power of features [10] and [11].
            const float MaxRtDeviationMinutes = 5.0f;

            if (!float.IsNaN(result.RtDeviationMinutes))
            {
                // Primary path: use calibrated deviation from RecalibrateRtDeviations.
                // Already in observed minutes; accounts for iRT-to-minutes conversion.
                float deltaRt = MathF.Min(MathF.Abs(result.RtDeviationMinutes), MaxRtDeviationMinutes);
                fv.RtDeviationMinutes = deltaRt;
                fv.RtDeviationSquared = deltaRt * deltaRt;
            }
            else if (result.LibraryRetentionTime.HasValue && !float.IsNaN(result.ObservedApexRt))
            {
                // Fallback: direct subtraction — only valid when LibraryRetentionTime
                // is already in observed minutes (e.g. benchmark with DIA-NN ground truth).
                float deltaRt = MathF.Abs(
                    result.ObservedApexRt - (float)result.LibraryRetentionTime.Value);
                deltaRt = MathF.Min(deltaRt, MaxRtDeviationMinutes);
                fv.RtDeviationMinutes = deltaRt;
                fv.RtDeviationSquared = deltaRt * deltaRt;
            }
            else
            {
                fv.RtDeviationMinutes = MaxRtDeviationMinutes;
                fv.RtDeviationSquared = MaxRtDeviationMinutes * MaxRtDeviationMinutes;
            }

            // ── Mass accuracy features [13-15] ──────────────────────────────
            // Use absolute value of mean signed error — direction of calibration offset
            // isn't discriminative, magnitude is.
            fv.MeanMassErrorPpm = float.IsNaN(result.MeanMassErrorPpm)
                ? 0f : MathF.Abs(result.MeanMassErrorPpm);
            fv.MassErrorStdPpm = SafeScore(result.MassErrorStdPpm);
            fv.MaxAbsMassErrorPpm = SafeScore(result.MaxAbsMassErrorPpm);

            // ── Best-fragment reference curve [16-19] ───────────────────────
            fv.BestFragCorrelationSum = SafeScore(result.BestFragCorrelationSum);
            fv.MedianFragRefCorr = SafeScore(result.MedianFragRefCorr);
            fv.MinFragRefCorr = float.IsNaN(result.MinFragRefCorr)
                ? -1f : result.MinFragRefCorr;
            fv.StdFragRefCorr = SafeScore(result.StdFragRefCorr);

            // ── Signal ratio deviation [20-22] (lower = better for targets) ─
            fv.MeanSignalRatioDeviation = SafeScore(result.MeanSignalRatioDev);
            fv.MaxSignalRatioDeviation = SafeScore(result.MaxSignalRatioDev);
            fv.StdSignalRatioDeviation = SafeScore(result.StdSignalRatioDev);

            // ── Smoothed correlations [23-24] ────────────────────────────────
            fv.SmoothedMeanFragCorr = SafeScore(result.SmoothedMeanFragCorr);
            fv.SmoothedMinFragCorr = float.IsNaN(result.SmoothedMinFragCorr)
                ? -1f : result.SmoothedMinFragCorr;

            // ── Signal-to-noise [25] ─────────────────────────────────────────
            fv.Log2SignalToNoise = SafeScore(result.Log2SignalToNoise);

            // ── Migrated features [26-28] ────────────────────────────────────
            fv.BestFragWeightedCosine = SafeScore(result.BestFragWeightedCosine);
            fv.BoundarySignalRatio = SafeScore(result.BoundarySignalRatio);
            fv.ApexToMeanRatio = SafeScore(result.ApexToMeanRatio);

            // ── MS1 features [29-32] ─────────────────────────────────────────
            // If index is null or has no MS1 scans, all four remain NaN on the result.
            // The classifier handles NaN via mean imputation.
            if (index != null && index.Ms1ScanCount > 0)
            {
                DiaMs1FeatureComputer.ComputeMs1Features(
                    result, index, bestFragXic, bestFragXicRts);
            }

            // NaN is passed through; WriteTo() writes as-is.
            fv.PrecursorXicApexIntensity = result.PrecursorXicApexIntensity;
            fv.IsotopePatternScore = result.IsotopePatternScore;
            fv.Ms1Ms2Correlation = result.Ms1Ms2Correlation;
            fv.PrecursorElutionScore = result.PrecursorElutionScore;

            // ── ChimericScore [33] ───────────────────────────────────────────
            // Computed during assembly by DiaLibraryQueryGenerator.ComputeChimericScores().
            // NaN if only one precursor exists in the window.
            fv.ChimericScore = result.ChimericScore;

            // ── LibraryCoverageFraction [34] ─────────────────────────────────
            // Intensity-weighted fraction of library fragments detected.
            fv.LibraryCoverageFraction = result.LibraryCoverageFraction;

            // ── Metadata ─────────────────────────────────────────────────────
            fv.ChargeState = result.ChargeState;
            fv.PrecursorMz = (float)result.PrecursorMz;
            fv.IsDecoy = result.IsDecoy;
            fv.PrecursorIndex = precursorIndex;

            return fv;
        }

        /// <summary>Replaces NaN scores with 0 for classifier safety.</summary>
        private static float SafeScore(float score) =>
            float.IsNaN(score) ? 0f : score;

        // ── Best-fragment XIC extraction ──────────────────────────────────────
        // Used by Phase14BenchmarkRunner and PostDiaSearchAnalysisTask to wire
        // Ms1Ms2Correlation (feature slot [31]).

        /// <summary>
        /// Extracts per-scan intensities for a single fragment m/z within an RT window.
        /// Returns empty arrays when the window has no scans or the fragment is out of range.
        /// </summary>
        public static void ExtractBestFragmentXic(
            DiaScanIndex index,
            float fragmentMz,
            int windowId,
            float rtMin,
            float rtMax,
            float ppmTolerance,
            out float[] intensities,
            out float[] rts)
        {
            if (fragmentMz <= 0f || rtMin > rtMax ||
                !index.TryGetScanRangeForWindow(windowId, out int winStart, out int winCount) ||
                winCount == 0)
            {
                intensities = Array.Empty<float>();
                rts = Array.Empty<float>();
                return;
            }

            float daltonTol = fragmentMz * ppmTolerance * 1e-6f;
            float mzLo = fragmentMz - daltonTol;
            float mzHi = fragmentMz + daltonTol;

            // First pass: count scans in RT window
            int windowCount = 0;
            for (int i = winStart; i < winStart + winCount; i++)
            {
                float rt = index.GetScanRt(i);
                if (rt < rtMin) continue;
                if (rt > rtMax) break;
                windowCount++;
            }

            if (windowCount == 0)
            {
                intensities = Array.Empty<float>();
                rts = Array.Empty<float>();
                return;
            }

            var rtArr = new float[windowCount];
            var intArr = new float[windowCount];
            int written = 0;

            for (int i = winStart; i < winStart + winCount && written < windowCount; i++)
            {
                float scanRt = index.GetScanRt(i);
                if (scanRt < rtMin) continue;
                if (scanRt > rtMax) break;

                rtArr[written] = scanRt;
                intArr[written] = SumInMzWindow(
                    index.GetScanMzSpan(i),
                    index.GetScanIntensitySpan(i),
                    mzLo, mzHi);
                written++;
            }

            rts = rtArr;
            intensities = intArr;
        }

        private static float SumInMzWindow(
            ReadOnlySpan<float> mzs,
            ReadOnlySpan<float> intensities,
            float mzLo,
            float mzHi)
        {
            if (mzs.IsEmpty) return 0f;

            // Binary search for the left edge
            int lo = 0, hi = mzs.Length;
            while (lo < hi)
            {
                int mid = lo + ((hi - lo) >> 1);
                if (mzs[mid] < mzLo) lo = mid + 1;
                else hi = mid;
            }

            float sum = 0f;
            for (int i = lo; i < mzs.Length && mzs[i] <= mzHi; i++)
                sum += intensities[i];
            return sum;
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    //  DiaMs1FeatureComputer
    //  Computes MS1-based features [29-32] for a single DiaSearchResult.
    //
    //  Phase 16B, Prompt 5: Initial implementation (averagine dot product).
    //  Phase 16B, Prompt 8: Fix 2 (median normalization) + Fix 3 (noise-gated
    //                        dot product) applied — both features still inverted.
    //  Phase 16B, Prompt 9: Revised fixes:
    //    [29] median → lower-quartile (Q25) background reference.
    //    [30] averagine dot product → simple M+1/M0 ratio at apex.
    //
    //  Called from DiaFeatureExtractor.ComputeFeatures() when a DiaScanIndex
    //  with MS1 scans is available. All four features default to NaN on
    //  DiaSearchResult; if MS1 extraction finds insufficient signal the NaN
    //  is preserved and the classifier applies mean imputation.
    // ═══════════════════════════════════════════════════════════════════════
    public static class DiaMs1FeatureComputer
    {
        // MS1 extraction tolerance — wider than MS2 to account for lower
        // resolving power of survey scans and minor calibration drift.
        private const float DefaultMs1PpmTolerance = 20f;

        // Minimum MS1 scan points required for meaningful feature computation.
        private const int MinMs1Points = 3;

        /// <summary>
        /// Computes all four MS1 features and stores them on <paramref name="result"/>.
        /// If the index has no MS1 scans, all four features remain at their NaN defaults.
        /// </summary>
        /// <param name="result">
        /// Scored DiaSearchResult. Must have RtWindowStart/End, ChargeState,
        /// PrecursorMz, and BestFragIndex populated by the MS2 scoring step.
        /// </param>
        /// <param name="index">DiaScanIndex containing MS1 scan data.</param>
        /// <param name="bestFragXic">
        /// Intensities of the best MS2 fragment XIC over the RT window.
        /// Pass default to skip Ms1Ms2Correlation.
        /// </param>
        /// <param name="bestFragXicRts">
        /// Retention times corresponding to <paramref name="bestFragXic"/>.
        /// Pass default to skip Ms1Ms2Correlation.
        /// </param>
        /// <param name="ppmTolerance">MS1 extraction tolerance in ppm (default 20).</param>
        public static void ComputeMs1Features(
            DiaSearchResult result,
            DiaScanIndex index,
            ReadOnlySpan<float> bestFragXic = default,
            ReadOnlySpan<float> bestFragXicRts = default,
            float ppmTolerance = DefaultMs1PpmTolerance)
        {
            if (result == null) throw new ArgumentNullException(nameof(result));
            if (index == null) throw new ArgumentNullException(nameof(index));

            if (index.Ms1ScanCount == 0) return; // all four remain NaN

            if (ppmTolerance <= 0f) ppmTolerance = DefaultMs1PpmTolerance;

            float rtMin = result.RtWindowStart;
            float rtMax = result.RtWindowEnd;
            float precursorMz = (float)result.PrecursorMz;
            int chargeState = Math.Max(1, result.ChargeState);

            // ← ADD THESE TWO LINES:
            if (result.PrecursorMz == 0)  // only print for the very first precursor
                Console.WriteLine($"  DEBUG MS1[0]: rtMin={rtMin:F4} rtMax={rtMax:F4} mz={precursorMz:F4} z={chargeState} ms1Count={index.Ms1ScanCount}");


            // Extract M0, M+1, M+2 isotope XICs in one pass
            Ms1XicExtractor.ExtractIsotopeXics(
                index, precursorMz, chargeState,
                rtMin, rtMax, ppmTolerance,
                out float[] xicRts,
                out float[] m0Int,
                out float[] m1Int,
                out float[] m2Int);

            if (xicRts.Length < MinMs1Points) return;

            // ── [29] PrecursorXicApexIntensity ───────────────────────────────
            // log2(apexM0 / Q25(non-zero M0 XIC))
            // Q25 = genuine background floor; not inflated by the elution peak body.
            // Targets: apex ≫ Q25 (sharp peak) → high ratio.
            // Decoys:  no real peak, flat XIC → apex ≈ Q25 → low ratio.
            float apexM0 = FindMax(m0Int);
            if (apexM0 > 0f)
            {
                float q25M0 = LowerQuartileNonZero(m0Int);
                if (q25M0 > 0f)
                    result.PrecursorXicApexIntensity = MathF.Log2(apexM0 / q25M0);
            }

            // ── [30] IsotopePatternScore ──────────────────────────────────────
            // M+1/M0 intensity ratio at the XIC apex, clamped to [0, 1].
            // Noise gate: M+1 must exceed 5% of M0 to avoid noise-spike bias.
            int apexIdx = FindMaxIndex(m0Int);
            float obsM0 = m0Int[apexIdx];
            float obsM1 = m1Int[apexIdx];

            if (obsM0 > 0f && obsM1 > 0.05f * obsM0)
                result.IsotopePatternScore = Math.Clamp(obsM1 / obsM0, 0f, 1f);

            // ── [31] Ms1Ms2Correlation ────────────────────────────────────────
            // Pearson r between M0 XIC and best fragment XIC.
            if (!bestFragXic.IsEmpty && !bestFragXicRts.IsEmpty &&
                bestFragXic.Length == bestFragXicRts.Length)
            {
                float r = ComputeMs1Ms2Correlation(
                    xicRts, m0Int, bestFragXicRts, bestFragXic);
                if (!float.IsNaN(r))
                    result.Ms1Ms2Correlation = r;
            }

            // ── [32] PrecursorElutionScore ────────────────────────────────────
            // Gaussian fit quality of the M0 XIC.
            float gaussScore = ComputeGaussianFitScore(xicRts, m0Int);
            if (!float.IsNaN(gaussScore))
                result.PrecursorElutionScore = gaussScore;
        }

        // ── Private helpers ───────────────────────────────────────────────────

        /// <summary>
        /// Pearson correlation between the M0 precursor XIC and the best MS2 fragment XIC.
        /// XICs are sampled at different time points so we align via nearest-neighbor lookup.
        /// Returns NaN if fewer than MinMs1Points co-detected points exist after alignment.
        /// </summary>
        private static float ComputeMs1Ms2Correlation(
            ReadOnlySpan<float> ms1Rts,
            ReadOnlySpan<float> ms1Intensities,
            ReadOnlySpan<float> ms2Rts,
            ReadOnlySpan<float> ms2Intensities)
        {
            const float maxRtGapMin = 0.05f; // 3 s — handles typical ~2 s MS1 cycle

            int pairCount = 0;
            float sumA = 0f, sumB = 0f, sumAB = 0f, sumA2 = 0f, sumB2 = 0f;
            int ms1Cursor = 0;

            for (int j = 0; j < ms2Rts.Length; j++)
            {
                float ms2Rt = ms2Rts[j];
                float ms2Int = ms2Intensities[j];
                if (ms2Int <= 0f) continue;

                // Advance cursor to nearest MS1 scan (both arrays are RT-sorted)
                while (ms1Cursor + 1 < ms1Rts.Length &&
                       MathF.Abs(ms1Rts[ms1Cursor + 1] - ms2Rt) <
                       MathF.Abs(ms1Rts[ms1Cursor] - ms2Rt))
                    ms1Cursor++;

                if (MathF.Abs(ms1Rts[ms1Cursor] - ms2Rt) > maxRtGapMin) continue;

                float ms1Int = ms1Intensities[ms1Cursor];
                if (ms1Int <= 0f) continue;

                sumA += ms1Int;
                sumB += ms2Int;
                sumAB += ms1Int * ms2Int;
                sumA2 += ms1Int * ms1Int;
                sumB2 += ms2Int * ms2Int;
                pairCount++;
            }

            if (pairCount < MinMs1Points) return float.NaN;

            float denom = (pairCount * sumA2 - sumA * sumA) *
                          (pairCount * sumB2 - sumB * sumB);
            if (denom <= 0f) return float.NaN;

            return Math.Clamp(
                (pairCount * sumAB - sumA * sumB) / MathF.Sqrt(denom),
                -1f, 1f);
        }

        /// <summary>
        /// Gaussian fit score: cosine similarity between the observed XIC and a
        /// Gaussian predicted from the XIC's own weighted moments (μ, σ).
        /// Range [0, 1]; returns NaN if fewer than MinMs1Points or σ = 0.
        /// </summary>
        private static float ComputeGaussianFitScore(
            ReadOnlySpan<float> rts,
            ReadOnlySpan<float> intensities)
        {
            if (rts.Length < MinMs1Points) return float.NaN;

            // Step 1: weighted moments → μ, σ
            float sumW = 0f, sumWRt = 0f, sumWRt2 = 0f;
            for (int i = 0; i < rts.Length; i++)
            {
                float w = intensities[i];
                if (w <= 0f) continue;
                sumW += w;
                sumWRt += w * rts[i];
                sumWRt2 += w * rts[i] * rts[i];
            }

            if (sumW <= 0f) return float.NaN;

            float mu = sumWRt / sumW;
            float variance = sumWRt2 / sumW - mu * mu;
            if (variance <= 0f) return float.NaN;

            float sigma = MathF.Sqrt(variance);
            float apexObs = FindMax(intensities);
            if (apexObs <= 0f) return float.NaN;

            // Step 2: predicted Gaussian scaled to observed apex
            float[] predBuf = rts.Length > 256
                ? ArrayPool<float>.Shared.Rent(rts.Length)
                : null;

            Span<float> pred = rts.Length <= 256
                ? stackalloc float[rts.Length]
                : predBuf.AsSpan(0, rts.Length);

            float invSigma2 = 1f / (2f * sigma * sigma);
            for (int i = 0; i < rts.Length; i++)
            {
                float dt = rts[i] - mu;
                pred[i] = apexObs * MathF.Exp(-dt * dt * invSigma2);
            }

            // Step 3: cosine similarity
            float dot = 0f, normObs = 0f, normPred = 0f;
            for (int i = 0; i < rts.Length; i++)
            {
                dot += intensities[i] * pred[i];
                normObs += intensities[i] * intensities[i];
                normPred += pred[i] * pred[i];
            }

            if (predBuf != null) ArrayPool<float>.Shared.Return(predBuf);

            if (normObs <= 0f || normPred <= 0f) return float.NaN;

            return Math.Clamp(
                dot / (MathF.Sqrt(normObs) * MathF.Sqrt(normPred)),
                0f, 1f);
        }

        /// <summary>
        /// Returns the 25th percentile of all strictly-positive values in the span.
        /// Returns 0f if there are no positive values.
        /// Uses stackalloc for ≤ 256 elements, ArrayPool for larger arrays.
        /// </summary>
        private static float LowerQuartileNonZero(ReadOnlySpan<float> values)
        {
            int count = 0;
            for (int i = 0; i < values.Length; i++)
                if (values[i] > 0f) count++;

            if (count == 0) return 0f;

            float[] rentedBuf = count > 256 ? ArrayPool<float>.Shared.Rent(count) : null;
            Span<float> buf = count <= 256
                ? stackalloc float[count]
                : rentedBuf.AsSpan(0, count);

            try
            {
                int j = 0;
                for (int i = 0; i < values.Length; i++)
                    if (values[i] > 0f) buf[j++] = values[i];

                buf.Sort();

                // Nearest-rank Q25: index = floor(0.25 * count), clamped
                int q25Idx = Math.Max(0, (int)MathF.Floor(0.25f * count) - 1);
                return buf[q25Idx];
            }
            finally
            {
                if (rentedBuf != null) ArrayPool<float>.Shared.Return(rentedBuf);
            }
        }

        private static float FindMax(ReadOnlySpan<float> values)
        {
            float max = 0f;
            for (int i = 0; i < values.Length; i++)
                if (values[i] > max) max = values[i];
            return max;
        }

        private static int FindMaxIndex(ReadOnlySpan<float> values)
        {
            int idx = 0;
            float max = float.MinValue;
            for (int i = 0; i < values.Length; i++)
            {
                if (values[i] > max) { max = values[i]; idx = i; }
            }
            return idx;
        }
    }
}