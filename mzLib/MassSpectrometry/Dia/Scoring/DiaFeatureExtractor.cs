// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Buffers;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Feature vector for a single DIA precursor identification.
    /// 
    /// Phase 16A, Prompt 1: All 3 previously-unmigrated features now active.
    ///   ClassifierFeatureCount = 29 (26 original + 3 migrated).
    /// Phase 16B, Prompt 5-6: Four MS1 features added; ClassifierFeatureCount = 33.
    /// Phase 16B, Prompt 8: MS1 fixes applied (median norm + noise-gated dot product).
    ///   Both [29] and [30] remained inverted after Prompt 8.
    /// Phase 16B, Prompt 9: MS1 fixes revised:
    ///   [29] median → Q25 background reference (LowerQuartileNonZero).
    ///   [30] averagine dot product → simple M+1/M0 ratio.
    /// Phase 23 Prompt 10: PrecursorXicApexIntensity gating was attempted and reverted.
    ///   See DiaMs1FeatureComputer class comment for details.
    /// Prompt 3 (rebuild): CoElutionStd [35] and CandidateScoreGap [36] added.
    ///   ClassifierFeatureCount = 35 → 37.
    /// MS1 Interference Resolution phase: Ms1ApexConfirmationScore [37] added.
    ///   ClassifierFeatureCount = 37 → 38.
    /// 
    /// Migrated features (computed in AssembleResultsWithTemporalScoring but previously
    /// not forwarded to the classifier):
    ///   [26] BestFragWeightedCosine  — weighted cosine using best-fragment correlations as weights
    ///   [27] BoundarySignalRatio     — boundary TIC / apex TIC (low = sharp peak)
    ///   [28] ApexToMeanRatio         — apex TIC / mean TIC in peak range (high = prominent peak)
    /// 
    /// MS1 features (wired into vector in Prompt 6, revised in Prompt 9/10):
    ///   [29] PrecursorXicApexIntensity — log2(apex M0 / Q25 M0) — peak prominence above background
    ///   [30] IsotopePatternScore       — M+1/M0 ratio at apex (noise-gated at 5% of M0)
    ///   [31] Ms1Ms2Correlation         — Pearson r between M0 XIC and best fragment XIC
    ///   [32] PrecursorElutionScore     — Gaussian fit quality of precursor XIC
    /// 
    /// Peak selection quality features (Prompt 3):
    ///   [35] CoElutionStd      — std dev of per-fragment apex positions (low = tight co-elution)
    ///   [36] CandidateScoreGap — best minus second-best candidate score (high = unambiguous peak)
    /// 
    /// Current layout: 37 classifier features.
    /// </summary>
    public struct DiaFeatureVector
    {
        // -- Primary scores (3) -------------------------------------------
        public float ApexScore;             // [0]
        public float TemporalScore;         // [1]
        public float SpectralAngle;         // [2]

        // -- Peak-region correlations (2) ---------------------------------
        public float PeakMeanFragCorr;      // [3]
        public float PeakMinFragCorr;       // [4]

        // -- Peak shape (2) -----------------------------------------------
        public float PeakWidth;             // [5]
        public float CandidateCount;        // [6]

        // -- Signal features (2) ------------------------------------------
        public float LogTotalIntensity;     // [7]
        public float IntensityCV;           // [8]

        // -- Fragment evidence (1) ----------------------------------------
        public float FragmentDetectionRate; // [9]

        // -- Retention time (2) -------------------------------------------
        public float RtDeviationMinutes;    // [10]
        public float RtDeviationSquared;    // [11]

        // -- Temporal evidence (1) ----------------------------------------
        public int TimePointsUsed;          // [12]

        // -- Mass accuracy (3) --------------------------------------------
        public float MeanMassErrorPpm;      // [13]  (absolute value of mean signed error)
        public float MassErrorStdPpm;       // [14]
        public float MaxAbsMassErrorPpm;    // [15]

        // -- Best-fragment reference curve (4) ----------------------------
        public float BestFragCorrelationSum;// [16]
        public float MedianFragRefCorr;     // [17]
        public float MinFragRefCorr;        // [18]
        public float StdFragRefCorr;        // [19]

        // -- Signal ratio deviation (3) -----------------------------------
        public float MeanSignalRatioDeviation;  // [20]
        public float MaxSignalRatioDeviation;   // [21]
        public float StdSignalRatioDeviation;   // [22]

        // -- Smoothed correlations (2) ------------------------------------
        public float SmoothedMeanFragCorr;  // [23]
        public float SmoothedMinFragCorr;   // [24]

        // -- Signal-to-noise (1) ------------------------------------------
        public float Log2SignalToNoise;     // [25]

        // -- Migrated features (3) -- Phase 16A ---------------------------
        public float BestFragWeightedCosine;   // [26]
        public float BoundarySignalRatio;      // [27]
        public float ApexToMeanRatio;          // [28]

        // -- MS1 features (4) -- Phase 16B, Prompt 6 ---------------------
        public float PrecursorXicApexIntensity; // [29]
        public float IsotopePatternScore;        // [30]
        public float Ms1Ms2Correlation;          // [31]
        public float PrecursorElutionScore;      // [32]

        // -- Interference / chimeric features (1) -- Phase 19, Priority 2
        public float ChimericScore;              // [33]

        // -- Coverage features (1) -- Phase 19, Priority 5
        public float LibraryCoverageFraction;    // [34]

        // -- Peak selection quality (2) -- Prompt 3 -----------------------
        /// <summary>
        /// Std dev of per-fragment apex positions within the detected peak window.
        /// Low value = fragments co-elute tightly (good). High = scattered (interference).
        /// Populated via DiaSearchResult.CoElutionStd.
        /// NaN default → sentinel 99f in WriteTo.
        /// </summary>
        public float CoElutionStd;               // [35]

        /// <summary>
        /// SelectionScore of the best peak group candidate minus the second-best.
        /// Large gap = unambiguous peak. Small gap = competing peaks present.
        /// NaN default → sentinel 0f in WriteTo.
        /// </summary>
        public float CandidateScoreGap;          // [36]

        // -- MS1 apex confirmation (1) -- MS1 Interference Resolution phase
        /// <summary>
        /// Raw MS1 apex intensity ratio for the selected peak group.
        /// Ratio of precursor MS1 intensity at the selected apex RT to the maximum precursor
        /// MS1 intensity in the search window. Range [0, 1].
        ///
        ///   1.0 = precursor has strong MS1 signal at the apex → strong target confirmation
        ///   0.0 = precursor has no MS1 signal at the apex → possible interference
        ///
        /// Default 1.0f (neutral) when MS1 is unavailable or window ≤ 1.0 min.
        /// Targets should have systematically higher values than decoys.
        /// </summary>
        public float Ms1ApexConfirmationScore;   // [37]

        // -- Metadata (not classifier features) ---------------------------
        public bool IsDecoy;
        public int PrecursorIndex;
        public int ChargeState;
        public float PrecursorMz;
        public int FragmentsDetected;
        public int FragmentsQueried;

        /// <summary>
        /// Number of features used by the classifier.
        /// Prompt 3: 35 → 37 (added CoElutionStd [35] and CandidateScoreGap [36]).
        /// MS1 Interference Resolution phase: 37 → 38 (added Ms1ApexConfirmationScore [37]).
        /// </summary>
        public const int ClassifierFeatureCount = 38;

        /// <summary>Index of ApexScore in the feature vector.</summary>
        public const int InteractionFeatureIndexA = 0;  // ApexScore

        /// <summary>Index of PeakMeanFragCorr in the feature vector.</summary>
        public const int InteractionFeatureIndexB = 3;  // PeakMeanFragCorr

        /// <summary>
        /// Writes classifier features into a float span.
        /// Order must be consistent with weight vectors.
        /// </summary>
        public readonly void WriteTo(Span<float> features)
        {
            if (features.Length < ClassifierFeatureCount)
                throw new ArgumentException($"Span must have at least {ClassifierFeatureCount} elements");

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
            features[26] = BestFragWeightedCosine;
            features[27] = BoundarySignalRatio;
            features[28] = ApexToMeanRatio;
            features[29] = PrecursorXicApexIntensity;
            features[30] = IsotopePatternScore;
            features[31] = Ms1Ms2Correlation;
            features[32] = PrecursorElutionScore;
            features[33] = ChimericScore;
            features[34] = LibraryCoverageFraction;
            // Peak selection quality [35-36]
            features[35] = float.IsNaN(CoElutionStd) ? 99f : CoElutionStd;
            features[36] = float.IsNaN(CandidateScoreGap) ? 0f : CandidateScoreGap;
            // MS1 apex confirmation [37]
            features[37] = float.IsNaN(Ms1ApexConfirmationScore) ? 1.0f : Ms1ApexConfirmationScore;
        }

        /// <summary>
        /// Feature names for TSV/diagnostic output. Order matches WriteTo().
        /// </summary>
        public static readonly string[] FeatureNames = new[]
        {
            "ApexScore", "TemporalScore", "SpectralAngle",                        // [0-2]
            "PeakMeanFragCorr", "PeakMinFragCorr",                                // [3-4]
            "PeakWidth", "CandidateCount",                                        // [5-6]
            "LogTotalIntensity", "IntensityCV",                                   // [7-8]
            "FragDetRate",                                                         // [9]
            "RtDeviationMinutes", "RtDeviationSquared",                           // [10-11]
            "TimePointsUsed",                                                     // [12]
            "MeanMassErrorPpm", "MassErrorStdPpm", "MaxAbsMassErrorPpm",          // [13-15]
            "BestFragCorrSum", "MedianFragRefCorr", "MinFragRefCorr", "StdFragRefCorr", // [16-19]
            "MeanSigRatioDev", "MaxSigRatioDev", "StdSigRatioDev",                // [20-22]
            "SmoothedMeanFragCorr", "SmoothedMinFragCorr",                        // [23-24]
            "Log2SNR",                                                            // [25]
            "BestFragWeightedCosine", "BoundarySignalRatio", "ApexToMeanRatio",   // [26-28]
            "PrecursorXicApexIntensity", "IsotopePatternScore",                   // [29-30]
            "Ms1Ms2Correlation", "PrecursorElutionScore",                         // [31-32]
            "ChimericScore",                                                      // [33]
            "LibraryCoverageFraction",                                            // [34]
            "CoElutionStd",                                                       // [35]
            "CandidateScoreGap",                                                  // [36]
            "Ms1ApexConfirmationScore",                                           // [37]
        };
    }

    /// <summary>
    /// Computes feature vectors from DiaSearchResult objects.
    /// Thread-safe: uses only local state + ArrayPool.
    /// </summary>
    public static class DiaFeatureExtractor
    {
        public static DiaFeatureVector ComputeFeatures(
            DiaSearchResult result,
            int precursorIndex,
            DiaScanIndex index = null,
            ReadOnlySpan<float> bestFragXic = default,
            ReadOnlySpan<float> bestFragXicRts = default)
        {
            var fv = new DiaFeatureVector();

            // -- Primary scores [0-2] ------------------------------------
            fv.ApexScore = SafeScore(result.ApexScore);
            fv.TemporalScore = SafeScore(result.TemporalScore);
            float _sa = result.SpectralAngle;
            fv.SpectralAngle = SafeScore(float.IsNaN(_sa) ? result.SpectralAngleScore : _sa);

            // -- Peak-region correlations [3-4] ---------------------------
            if (result.DetectedPeakGroup.HasValue && result.DetectedPeakGroup.Value.IsValid)
            {
                fv.PeakMeanFragCorr = SafeScore(result.PeakMeanFragCorr);
                fv.PeakMinFragCorr = float.IsNaN(result.PeakMinFragCorr)
                    ? -1f : result.PeakMinFragCorr;
                fv.PeakWidth = result.DetectedPeakGroup.Value.PeakWidthMinutes;
            }
            else
            {
                fv.PeakMeanFragCorr = SafeScore(result.MeanFragCorr);
                fv.PeakMinFragCorr = float.IsNaN(result.MinFragCorr)
                    ? -1f : result.MinFragCorr;
                fv.PeakWidth = 0f;
            }

            // -- Peak shape [5-6] ----------------------------------------
            fv.CandidateCount = SafeScore(result.CandidateCount);

            // -- Fragment evidence [9] ------------------------------------
            fv.FragmentDetectionRate = result.FragmentDetectionRate;
            fv.FragmentsDetected = result.FragmentsDetected;
            fv.FragmentsQueried = result.FragmentsQueried;

            // -- Temporal evidence [12] -----------------------------------
            fv.TimePointsUsed = result.TimePointsUsed;

            // -- Intensity features [7-8] ---------------------------------
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

            // -- Retention time features [10-11] --------------------------
            const float MaxRtDeviationMinutes = 5.0f;

            if (!float.IsNaN(result.RtDeviationMinutes))
            {
                float deltaRt = MathF.Min(MathF.Abs(result.RtDeviationMinutes), MaxRtDeviationMinutes);
                fv.RtDeviationMinutes = deltaRt;
                fv.RtDeviationSquared = deltaRt * deltaRt;
            }
            else
            {
                fv.RtDeviationMinutes = MaxRtDeviationMinutes;
                fv.RtDeviationSquared = MaxRtDeviationMinutes * MaxRtDeviationMinutes;
            }

            // -- Mass accuracy features [13-15] ---------------------------
            fv.MeanMassErrorPpm = float.IsNaN(result.MeanMassErrorPpm)
                ? 0f : MathF.Abs(result.MeanMassErrorPpm);
            fv.MassErrorStdPpm = SafeScore(result.MassErrorStdPpm);
            fv.MaxAbsMassErrorPpm = SafeScore(result.MaxAbsMassErrorPpm);

            // -- Best-fragment reference curve [16-19] --------------------
            fv.BestFragCorrelationSum = SafeScore(result.BestFragCorrelationSum);
            fv.MedianFragRefCorr = SafeScore(result.MedianFragRefCorr);
            fv.MinFragRefCorr = float.IsNaN(result.MinFragRefCorr)
                ? -1f : result.MinFragRefCorr;
            fv.StdFragRefCorr = SafeScore(result.StdFragRefCorr);

            // -- Signal ratio deviation [20-22] ---------------------------
            fv.MeanSignalRatioDeviation = SafeScore(result.MeanSignalRatioDev);
            fv.MaxSignalRatioDeviation = SafeScore(result.MaxSignalRatioDev);
            fv.StdSignalRatioDeviation = SafeScore(result.StdSignalRatioDev);

            // -- Smoothed correlations [23-24] ----------------------------
            fv.SmoothedMeanFragCorr = SafeScore(result.SmoothedMeanFragCorr);
            fv.SmoothedMinFragCorr = float.IsNaN(result.SmoothedMinFragCorr)
                ? -1f : result.SmoothedMinFragCorr;

            // -- Signal-to-noise [25] -------------------------------------
            fv.Log2SignalToNoise = SafeScore(result.Log2SignalToNoise);

            // -- Migrated features [26-28] --------------------------------
            fv.BestFragWeightedCosine = SafeScore(result.BestFragWeightedCosine);
            fv.BoundarySignalRatio = SafeScore(result.BoundarySignalRatio);
            fv.ApexToMeanRatio = SafeScore(result.ApexToMeanRatio);

            // -- MS1 features [29-32] -------------------------------------
            if (index != null && index.Ms1ScanCount > 0)
            {
                DiaMs1FeatureComputer.ComputeMs1Features(
                    result, index, bestFragXic, bestFragXicRts);
            }

            fv.PrecursorXicApexIntensity = result.PrecursorXicApexIntensity;
            fv.IsotopePatternScore = result.IsotopePatternScore;
            fv.Ms1Ms2Correlation = result.Ms1Ms2Correlation;
            fv.PrecursorElutionScore = result.PrecursorElutionScore;

            // -- ChimericScore [33] ---------------------------------------
            fv.ChimericScore = result.ChimericScore;

            // -- LibraryCoverageFraction [34] -----------------------------
            fv.LibraryCoverageFraction = result.LibraryCoverageFraction;

            // -- Peak selection quality [35-36] ---------------------------
            // CoElutionStd: initialized to 0f on DiaSearchResult; 0f = perfect co-elution (good).
            // CandidateScoreGap: initialized to 0f; 0f = no gap / single candidate.
            // WriteTo applies NaN sentinels (99f and 0f) for any NaN cases.
            fv.CoElutionStd = result.CoElutionStd;
            fv.CandidateScoreGap = result.CandidateScoreGap;

            // -- MS1 apex confirmation [37] --------------------------------
            // Default 1.0f (neutral) is the initialization on DiaSearchResult.
            // WriteTo uses 1.0f as the NaN sentinel (same neutral value).
            fv.Ms1ApexConfirmationScore = result.Ms1ApexConfirmationScore;

            // -- Metadata -------------------------------------------------
            fv.ChargeState = result.ChargeState;
            fv.PrecursorMz = (float)result.PrecursorMz;
            fv.IsDecoy = result.IsDecoy;
            fv.PrecursorIndex = precursorIndex;

            return fv;
        }

        /// <summary>Replaces NaN scores with 0 for classifier safety.</summary>
        private static float SafeScore(float score) =>
            float.IsNaN(score) ? 0f : score;

        /// <summary>
        /// Extracts per-scan intensities for a single fragment m/z within an RT window.
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
                intArr[written] = SumInMzWindow(index.GetScanMzSpan(i), index.GetScanIntensitySpan(i), mzLo, mzHi);
                written++;
            }

            rts = rtArr;
            intensities = intArr;
        }

        private static float SumInMzWindow(
            ReadOnlySpan<float> mzs, ReadOnlySpan<float> intensities,
            float mzLo, float mzHi)
        {
            if (mzs.IsEmpty) return 0f;
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

    /// <summary>
    /// Computes the four MS1 features for a single precursor result.
    ///   [29] PrecursorXicApexIntensity — log2(apex M0 / Q25 M0)
    ///   [30] IsotopePatternScore       — M+1/M0 ratio at apex
    ///   [31] Ms1Ms2Correlation         — Pearson r between M0 XIC and best fragment XIC
    ///   [32] PrecursorElutionScore     — Gaussian fit quality of precursor XIC
    /// </summary>
    public static class DiaMs1FeatureComputer
    {
        private const float DefaultMs1PpmTolerance = 20f;
        private const int MinMs1Points = 3;

        public static void ComputeMs1Features(
            DiaSearchResult result,
            DiaScanIndex index,
            ReadOnlySpan<float> bestFragXic,
            ReadOnlySpan<float> bestFragXicRts,
            float ppmTolerance = DefaultMs1PpmTolerance)
        {
            if (result == null) throw new ArgumentNullException(nameof(result));
            if (index == null) throw new ArgumentNullException(nameof(index));

            if (index.Ms1ScanCount == 0) return;

            if (ppmTolerance <= 0f) ppmTolerance = DefaultMs1PpmTolerance;

            float rtMin = result.RtWindowStart;
            float rtMax = result.RtWindowEnd;
            float precursorMz = (float)result.PrecursorMz;
            int chargeState = Math.Max(1, result.ChargeState);

            Ms1XicExtractor.ExtractIsotopeXics(
                index, precursorMz, chargeState,
                rtMin, rtMax, ppmTolerance,
                out float[] xicRts,
                out float[] m0Int,
                out float[] m1Int,
                out float[] m2Int);

            if (xicRts.Length < MinMs1Points)
                return;

            // ── [31] Ms1Ms2Correlation ────────────────────────────────────────
            if (!bestFragXic.IsEmpty && !bestFragXicRts.IsEmpty &&
                bestFragXic.Length == bestFragXicRts.Length)
            {
                float r = ComputeMs1Ms2Correlation(
                    xicRts, m0Int,
                    bestFragXicRts, bestFragXic);

                if (!float.IsNaN(r))
                    result.Ms1Ms2Correlation = r;
            }

            // ── [29] PrecursorXicApexIntensity ───────────────────────────────
            float apexM0 = FindMax(m0Int);
            if (apexM0 > 0f)
            {
                float q25M0 = LowerQuartileNonZero(m0Int);
                if (q25M0 > 0f)
                    result.PrecursorXicApexIntensity = MathF.Log2(apexM0 / q25M0);
            }

            // ── [30] IsotopePatternScore ──────────────────────────────────────
            int apexIdx = FindMaxIndex(m0Int);
            float obsM0 = m0Int[apexIdx];
            float obsM1 = m1Int[apexIdx];

            if (obsM0 > 0f && obsM1 > 0.05f * obsM0)
            {
                result.IsotopePatternScore = Math.Clamp(obsM1 / obsM0, 0f, 1f);
            }

            // ── [32] PrecursorElutionScore ────────────────────────────────────
            float gaussScore = ComputeGaussianFitScore(xicRts, m0Int);
            if (!float.IsNaN(gaussScore))
                result.PrecursorElutionScore = gaussScore;
        }

        // ── Feature sub-computations ──────────────────────────────────────────

        private static float ComputeMs1Ms2Correlation(
            ReadOnlySpan<float> ms1Rts,
            ReadOnlySpan<float> ms1Intensities,
            ReadOnlySpan<float> ms2Rts,
            ReadOnlySpan<float> ms2Intensities)
        {
            const float maxRtGapMin = 0.05f;

            int pairCount = 0;
            float sumA = 0f, sumB = 0f, sumAB = 0f, sumA2 = 0f, sumB2 = 0f;

            int ms1Cursor = 0;
            for (int j = 0; j < ms2Rts.Length; j++)
            {
                float ms2Rt = ms2Rts[j];
                float ms2Int = ms2Intensities[j];
                if (ms2Int <= 0f) continue;

                while (ms1Cursor + 1 < ms1Rts.Length &&
                       MathF.Abs(ms1Rts[ms1Cursor + 1] - ms2Rt) <
                       MathF.Abs(ms1Rts[ms1Cursor] - ms2Rt))
                {
                    ms1Cursor++;
                }

                float gap = MathF.Abs(ms1Rts[ms1Cursor] - ms2Rt);
                if (gap > maxRtGapMin) continue;

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

            float denom = (pairCount * sumA2 - sumA * sumA) * (pairCount * sumB2 - sumB * sumB);
            if (denom <= 0f) return float.NaN;

            float r = (pairCount * sumAB - sumA * sumB) / MathF.Sqrt(denom);
            return Math.Clamp(r, -1f, 1f);
        }

        private static float ComputeGaussianFitScore(
            ReadOnlySpan<float> rts,
            ReadOnlySpan<float> intensities)
        {
            if (rts.Length < MinMs1Points) return float.NaN;

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

            float[] predBuf = rts.Length <= 256
                ? null
                : ArrayPool<float>.Shared.Rent(rts.Length);

            Span<float> pred = rts.Length <= 256
                ? stackalloc float[rts.Length]
                : predBuf.AsSpan(0, rts.Length);

            float invSigma2 = 1f / (2f * sigma * sigma);
            for (int i = 0; i < rts.Length; i++)
            {
                float dt = rts[i] - mu;
                pred[i] = apexObs * MathF.Exp(-dt * dt * invSigma2);
            }

            float dot = 0f, normObs = 0f, normPred = 0f;
            for (int i = 0; i < rts.Length; i++)
            {
                dot += intensities[i] * pred[i];
                normObs += intensities[i] * intensities[i];
                normPred += pred[i] * pred[i];
            }

            if (predBuf != null)
                ArrayPool<float>.Shared.Return(predBuf);

            if (normObs <= 0f || normPred <= 0f) return float.NaN;

            float score = dot / (MathF.Sqrt(normObs) * MathF.Sqrt(normPred));
            return Math.Clamp(score, 0f, 1f);
        }

        // ── Span utilities ───────────────────────────────────────────────────

        private static float LowerQuartileNonZero(ReadOnlySpan<float> values)
        {
            int count = 0;
            for (int i = 0; i < values.Length; i++)
                if (values[i] > 0f) count++;

            if (count == 0) return 0f;

            float[] rentedBuf = count > 256 ? ArrayPool<float>.Shared.Rent(count) : null;
            Span<float> buf = count <= 256 ? stackalloc float[count] : rentedBuf.AsSpan(0, count);

            try
            {
                int j = 0;
                for (int i = 0; i < values.Length; i++)
                    if (values[i] > 0f) buf[j++] = values[i];

                buf.Sort();

                int q25Idx = Math.Max(0, (int)MathF.Floor(0.25f * count) - 1);
                return buf[q25Idx];
            }
            finally
            {
                if (rentedBuf != null)
                    ArrayPool<float>.Shared.Return(rentedBuf);
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
                if (values[i] > max)
                {
                    max = values[i];
                    idx = i;
                }
            }
            return idx;
        }
    }
}