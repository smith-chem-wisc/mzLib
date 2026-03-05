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
    /// 
    /// Migrated features (computed in AssembleResultsWithTemporalScoring but previously
    /// not forwarded to the classifier):
    ///   [26] BestFragWeightedCosine  — weighted cosine using best-fragment correlations as weights
    ///   [27] BoundarySignalRatio     — boundary TIC / apex TIC (low = sharp peak)
    ///   [28] ApexToMeanRatio         — apex TIC / mean TIC in peak range (high = prominent peak)
    /// 
    /// MS1 features (wired into vector in Prompt 6, revised in Prompt 9):
    ///   [29] PrecursorXicApexIntensity — log2(apex M0 / Q25 M0) — peak prominence above background
    ///   [30] IsotopePatternScore       — M+1/M0 ratio at apex (noise-gated at 5% of M0)
    ///   [31] Ms1Ms2Correlation         — Pearson r between M0 XIC and best fragment XIC
    ///   [32] PrecursorElutionScore     — Gaussian fit quality of precursor XIC
    /// 
    /// DROPPED (from Phase 12/early Phase 13):
    ///   - MeanFragCorr / MinFragCorr -- replaced by SmoothedMeanFragCorr/SmoothedMinFragCorr + best-fragment features
    ///   - PeakApexScore -- negative separation in Phase 12 analysis (hurts discrimination)
    ///   - PeakSymmetry -- zero separation (0.01)
    ///   - MedianXicDepth / XicDepthCV -- weak features, not in Prompt 8 plan
    /// 
    /// ADDED (Phase 13-14):
    ///   - CandidateCount (strong discriminator, separation 1.66)
    ///   - MeanMassErrorPpm (absolute), MassErrorStdPpm, MaxAbsMassErrorPpm (mass accuracy: 3)
    ///   - SmoothedMinFragCorr (was computed but not in vector)
    /// 
    /// Current layout: 33 classifier features (29 after Prompt 5 → 33 after Prompt 6).
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

        // -- Migrated features (3) -- Action Item 5 -----------------------
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
        // RtDeviationNormalized dropped: 100% NaN (PeakWidth=0 when no peak group detected).
        public float LibraryCoverageFraction;    // [34]

        // -- Metadata (not classifier features) ---------------------------
        public bool IsDecoy;
        public int PrecursorIndex;
        public int ChargeState;
        public float PrecursorMz;
        public int FragmentsDetected;
        public int FragmentsQueried;

        /// <summary>
        /// Number of features used by the classifier.
        /// Phase 16A, Prompt 1: 29 features (26 original + 3 migrated: [26-28]).
        /// Phase 16B, Prompt 6: 33 features (29 + 4 MS1: [29-32]).
        /// Phase 19: 35 features (33 + ChimericScore[33] + LibraryCoverageFraction[34]).
        ///   RtDeviationNormalized dropped: 100% NaN due to PeakWidth=0 when no peak group detected.
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
        /// Phase 16B layout (33 features):
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
        ///   [26]    Migrated: BestFragWeightedCosine
        ///   [27]    Migrated: BoundarySignalRatio
        ///   [28]    Migrated: ApexToMeanRatio
        ///   [29]    MS1: PrecursorXicApexIntensity
        ///   [30]    MS1: IsotopePatternScore
        ///   [31]    MS1: Ms1Ms2Correlation
        ///   [32]    MS1: PrecursorElutionScore
        ///   [33]    Interference: ChimericScore
        ///   [34]    Coverage: LibraryCoverageFraction
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
    /// provided, ComputeMs1Features() is called first to populate features [29-32] on
    /// the result, then those values are read back into the vector. When index is null
    /// or has no MS1 scans, features [29-32] remain NaN (handled by mean imputation
    /// in the classifier).
    /// </summary>
    public static class DiaFeatureExtractor
    {
        /// <summary>
        /// Computes a 33-feature vector from a scored DiaSearchResult.
        /// 
        /// MS1 features [29-32] require a DiaScanIndex built with MS1 scans.
        /// Pass null for <paramref name="index"/> to produce a 33-feature vector
        /// where [29-32] are NaN (e.g., for MS2-only files or unit tests).
        /// </summary>
        /// <param name="result">Scored result from the extraction pipeline.</param>
        /// <param name="precursorIndex">Zero-based index of this precursor in the results list.</param>
        /// <param name="index">
        /// Optional DiaScanIndex with MS1 scan data. When provided and Ms1ScanCount > 0,
        /// MS1 features are computed and stored on <paramref name="result"/> before being
        /// read into the vector.
        /// </param>
        /// <param name="bestFragXic">
        /// Optional: intensities of the best MS2 fragment XIC, used for Ms1Ms2Correlation.
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

            // -- Primary scores [0-2] ------------------------------------
            fv.ApexScore = SafeScore(result.ApexScore);
            fv.TemporalScore = SafeScore(result.TemporalScore);
            // Coalesce: temporal path sets SpectralAngle; simple assembler sets SpectralAngleScore
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
                // Fallback: use full-window correlations when no peak detected
                fv.PeakMeanFragCorr = SafeScore(result.MeanFragCorr);
                fv.PeakMinFragCorr = float.IsNaN(result.MinFragCorr)
                    ? -1f : result.MinFragCorr;
                fv.PeakWidth = 0f;
            }

            // -- Peak shape [5-6] ----------------------------------------
            // PeakWidth already set above
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

            if (result.LibraryRetentionTime.HasValue && !float.IsNaN(result.ObservedApexRt))
            {
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
            // Compute and store on result first, then read back into vector.
            // If index is null or has no MS1 scans, all four remain NaN on the result
            // and the vector fields are set to NaN (WriteTo preserves NaN for these slots).
            if (index != null && index.Ms1ScanCount > 0)
            {
                DiaMs1FeatureComputer.ComputeMs1Features(
                    result, index, bestFragXic, bestFragXicRts);
            }

            // NaN is passed through for these four: WriteTo() writes them as-is,
            // and the classifier applies mean imputation for NaN values.
            fv.PrecursorXicApexIntensity = result.PrecursorXicApexIntensity;
            fv.IsotopePatternScore = result.IsotopePatternScore;
            fv.Ms1Ms2Correlation = result.Ms1Ms2Correlation;
            fv.PrecursorElutionScore = result.PrecursorElutionScore;

            // -- ChimericScore [33] ---------------------------------------
            // Computed during assembly in DiaLibraryQueryGenerator and stored on result.
            // NaN if only one precursor exists in the window (no co-isolation possible).
            fv.ChimericScore = result.ChimericScore;

            // -- LibraryCoverageFraction [34] ------------------------------
            // Intensity-weighted fraction of library fragments detected.
            fv.LibraryCoverageFraction = result.LibraryCoverageFraction;

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

        // ── Best-fragment XIC extraction ──────────────────────────────────
        // Used by both Phase14BenchmarkRunner and PostDiaSearchAnalysisTask
        // to wire Ms1Ms2Correlation (feature slot [31]).

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
    /// Phase 16B, Prompt 5: Computes the four MS1 features for a single precursor result.
    /// Phase 16B, Prompt 8: Fix 2 (PrecursorXicApexIntensity median normalization) and
    ///                       Fix 3 (IsotopePatternScore noise-gated dot product) applied.
    ///                       Fix 1 (Ms1Ms2Correlation wiring) is in Phase14BenchmarkRunner.
    /// Phase 16B, Prompt 9: Fix 2 revised — median → lower-quartile (Q25) background reference.
    ///                       Fix 3 revised — dot product → simple M+1/M0 ratio feature.
    /// 
    /// Called after MS2 scoring (in AssembleResultsWithTemporalScoring or equivalent)
    /// when a DiaScanIndex with MS1 data is available. If the index has no MS1 scans,
    /// all four features remain NaN on the result (the classifier handles NaN via mean
    /// imputation; LDA is unaffected by a constant-NaN column).
    /// 
    /// The four features are stored on DiaSearchResult and will be wired into
    /// DiaFeatureVector [29-32] in Prompt 6.
    /// 
    /// Design notes:
    ///   - Uses the calibrated RT window (RtWindowStart/End) from DiaSearchResult, not a
    ///     hardcoded value, so extraction always matches the MS2 extraction window.
    ///   - All temporary arrays are rented from ArrayPool; no per-precursor heap allocation.
    ///   - BestFragXic is derived from DiaSearchResult.BestFragIndex and the pre-extracted
    ///     MS2 fragment XIC data stored on the result.
    /// 
    /// Prompt 8 changes (superseded by Prompt 9):
    ///   [29] PrecursorXicApexIntensity: log2(apexM0 / median(non-zero M0 XIC)).
    ///        Was still inverted post-Prompt 8 (target 1.511 < decoy 2.518).
    ///   [30] IsotopePatternScore: noise-gated averagine dot product.
    ///        Was still inverted post-Prompt 8 (decoy 0.459 > target 0.430).
    /// 
    /// Prompt 9 changes (this version):
    ///   [29] PrecursorXicApexIntensity: log2(apexM0 / Q25(non-zero M0 XIC)).
    ///        Lower quartile = genuine background floor; not inflated by the peak body.
    ///        Expected: inversion corrected, sep → ~0.5-1.0, +150-300 IDs.
    ///   [30] IsotopePatternScore: M+1/M0 ratio at apex, clamped [0,1], 5% noise gate.
    ///        Robust to averagine mismatch and co-eluter bias.
    ///        Expected: inversion corrected, sep → ~0.1-0.3, NaN rate < 20%, +50-150 IDs.
    /// </summary>
    public static class DiaMs1FeatureComputer
    {
        // Extraction tolerance in ppm for MS1 peak summing.
        // 20 ppm is wider than typical Orbitrap MS2 tolerance (10 ppm) to account for
        // the lower resolving power of survey scans and minor calibration drift.
        private const float DefaultMs1PpmTolerance = 20f;

        // Minimum number of MS1 scan points required for meaningful feature computation.
        private const int MinMs1Points = 3;

        /// <summary>
        /// Computes all four MS1 features and stores them on <paramref name="result"/>.
        /// 
        /// If the index has no MS1 scans, all four features remain at their NaN defaults
        /// and this method returns immediately.
        /// </summary>
        /// <param name="result">
        /// Scored DiaSearchResult. Must have RtWindowStart/End set (calibrated window),
        /// ChargeState, PrecursorMz, and BestFragIndex populated by the MS2 scoring step.
        /// </param>
        /// <param name="index">DiaScanIndex containing MS1 scan data.</param>
        /// <param name="bestFragXic">
        /// Intensities of the best MS2 fragment XIC over the RT window, aligned to the
        /// same time points as the MS1 scan RTs. Pass null or empty to skip Ms1Ms2Correlation.
        /// Typically derived from the flat intensity matrix column for BestFragIndex.
        /// </param>
        /// <param name="bestFragXicRts">
        /// Retention times corresponding to <paramref name="bestFragXic"/>. Must be parallel.
        /// Pass null or empty to skip Ms1Ms2Correlation.
        /// </param>
        /// <param name="ppmTolerance">
        /// MS1 extraction tolerance in ppm. Defaults to 20 ppm if ≤ 0.
        /// </param>
        public static void ComputeMs1Features(
            DiaSearchResult result,
            DiaScanIndex index,
            ReadOnlySpan<float> bestFragXic,
            ReadOnlySpan<float> bestFragXicRts,
            float ppmTolerance = DefaultMs1PpmTolerance)
        {
            if (result == null) throw new ArgumentNullException(nameof(result));
            if (index == null) throw new ArgumentNullException(nameof(index));

            // Early exit — all four features remain NaN
            if (index.Ms1ScanCount == 0) return;

            if (ppmTolerance <= 0f) ppmTolerance = DefaultMs1PpmTolerance;

            float rtMin = result.RtWindowStart;
            float rtMax = result.RtWindowEnd;
            float precursorMz = (float)result.PrecursorMz;
            int chargeState = Math.Max(1, result.ChargeState);

            // Extract isotope XICs (M0, M+1, M+2) in one pass
            Ms1XicExtractor.ExtractIsotopeXics(
                index, precursorMz, chargeState,
                rtMin, rtMax, ppmTolerance,
                out float[] xicRts,
                out float[] m0Int,
                out float[] m1Int,
                out float[] m2Int);

            if (xicRts.Length < MinMs1Points)
                return; // RT window too sparse for MS1 — all features remain NaN

            // ── [29] PrecursorXicApexIntensity ───────────────────────────────
            // Lower-quartile-normalized log2 apex: log2(apexM0 / Q25(non-zero M0 XIC)).
            //
            // Prompt 8 fix: median normalization was still inverted (target 1.511 < decoy 2.518).
            // Root cause: the ±0.30 min calibrated RT window is wide relative to peptide FWHM
            // (~5-10 s). For a true target, the apex scan is a large fraction of the sorted
            // intensity distribution, pulling the median upward and suppressing the apex/median
            // ratio. Decoys accumulate stray co-eluter signal spread across scans, giving a
            // flatter XIC where the median stays low and the ratio inflates.
            //
            // Prompt 9 fix: use the 25th percentile (lower quartile) as the background reference.
            // Q25 reflects the genuine low-signal background floor of the window, not the peak
            // body. For true targets: apex ≫ Q25 (sharp peak stands out from background) → 
            // high log2 ratio. For decoys: no real peak, XIC is flat noise → apex ≈ Q25 → 
            // low log2 ratio. This restores the correct target > decoy direction.
            //
            // Expected: sep flips positive (~0.5–1.0), +150–300 IDs.
            float apexM0 = FindMax(m0Int);
            if (apexM0 > 0f)
            {
                float q25M0 = LowerQuartileNonZero(m0Int);
                if (q25M0 > 0f)
                    result.PrecursorXicApexIntensity = MathF.Log2(apexM0 / q25M0);
                // else: all scans flat/zero — no useful relative signal, remains NaN
            }
            // else: no M0 signal at all, remains NaN

            // ── [30] IsotopePatternScore ──────────────────────────────────────
            // M+1/M0 intensity ratio at the XIC apex, clamped to [0, 1].
            //
            // Prompt 8 fix: noise-gated averagine dot product was still inverted
            // (decoy 0.459 > target 0.430, sep = 0.029). Root cause: the dot product
            // scores how well the observed envelope matches the averagine shape, but
            // co-eluting peptides systematically inflate M+1 for both targets and decoys.
            // Averagine is also a poor fit for this dataset's peptide mass range.
            //
            // Prompt 9 fix: replace with a simple M+1/M0 ratio at the apex scan.
            // Rationale:
            //   - A real peptide's M+1 is a predictable, non-zero fraction of M0.
            //     True target: apex is a real elution peak → M0 is large, M+1 is
            //     present at the correct isotopic fraction → ratio is stable and positive.
            //   - A decoy m/z is not a real peptide → M0 is noise; M+1 channel may be
            //     noise at a different level → ratio is low or erratic.
            //   - Ratio is insensitive to absolute intensity and gradient TIC variation.
            //   - No model-fitting: lower NaN rate expected vs dot product (30.8% → <15%).
            //   - Range [0, 1]: clamp prevents leverage from accidental high M+1 (co-eluter).
            //
            // Noise gate retained: M+1 must be > 5% of M0. Prevents noise spikes in
            // empty M+1 windows from producing a spurious near-zero ratio.
            // If below threshold: feature remains NaN (mean-imputed by classifier).
            //
            // Expected: inversion corrected, sep → ~0.1–0.3, NaN rate < 20%, +50–150 IDs.
            int apexIdx = FindMaxIndex(m0Int);
            float obsM0 = m0Int[apexIdx];
            float obsM1 = m1Int[apexIdx];

            if (obsM0 > 0f && obsM1 > 0.05f * obsM0)
            {
                // Clamp to [0, 1]: ratio > 1 is possible for a multiply-charged heavy
                // peptide or when a co-eluter strongly overlaps the M+1 channel.
                result.IsotopePatternScore = Math.Clamp(obsM1 / obsM0, 0f, 1f);
            }
            // else: M+1 below noise gate or no M0 signal → remains NaN

            // ── [31] Ms1Ms2Correlation ────────────────────────────────────────
            // Pearson r between M0 XIC and best fragment XIC over shared time points.
            if (!bestFragXic.IsEmpty && !bestFragXicRts.IsEmpty &&
                bestFragXic.Length == bestFragXicRts.Length)
            {
                float r = ComputeMs1Ms2Correlation(
                    xicRts, m0Int,
                    bestFragXicRts, bestFragXic);

                if (!float.IsNaN(r))
                    result.Ms1Ms2Correlation = r;
            }
            // else: remains NaN

            // ── [32] PrecursorElutionScore ────────────────────────────────────
            // Gaussian fit quality of the M0 XIC — how well it matches a Gaussian shape.
            float gaussScore = ComputeGaussianFitScore(xicRts, m0Int);
            if (!float.IsNaN(gaussScore))
                result.PrecursorElutionScore = gaussScore;
            // else: remains NaN
        }

        // ── Feature sub-computations ──────────────────────────────────────────

        /// <summary>
        /// Computes Pearson correlation between the M0 precursor XIC and the best MS2
        /// fragment XIC. The two XICs may be sampled at different time points (MS1 scans
        /// occur ~every 2s; MS2 windows every ~0.04 min per cycle), so we align them by
        /// interpolating the M0 intensity at each MS2 time point via nearest-neighbor lookup.
        /// 
        /// Returns NaN if fewer than MinMs1Points co-detected points exist after alignment.
        /// </summary>
        private static float ComputeMs1Ms2Correlation(
            ReadOnlySpan<float> ms1Rts,
            ReadOnlySpan<float> ms1Intensities,
            ReadOnlySpan<float> ms2Rts,
            ReadOnlySpan<float> ms2Intensities)
        {
            // Build aligned pairs: for each MS2 time point, find the nearest MS1 scan
            // within ±0.05 min (3 s). This handles typical ~2s MS1 cycle time.
            const float maxRtGapMin = 0.05f;

            int pairCount = 0;
            float sumA = 0f, sumB = 0f, sumAB = 0f, sumA2 = 0f, sumB2 = 0f;

            int ms1Cursor = 0; // advance forward as MS2 RTs increase
            for (int j = 0; j < ms2Rts.Length; j++)
            {
                float ms2Rt = ms2Rts[j];
                float ms2Int = ms2Intensities[j];
                if (ms2Int <= 0f) continue;

                // Advance ms1Cursor to the nearest MS1 scan to this MS2 RT
                // (MS1 RTs are sorted, so we can use a linear scan from the cursor)
                while (ms1Cursor + 1 < ms1Rts.Length &&
                       MathF.Abs(ms1Rts[ms1Cursor + 1] - ms2Rt) <
                       MathF.Abs(ms1Rts[ms1Cursor] - ms2Rt))
                {
                    ms1Cursor++;
                }

                float gap = MathF.Abs(ms1Rts[ms1Cursor] - ms2Rt);
                if (gap > maxRtGapMin) continue; // no MS1 scan close enough

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

        /// <summary>
        /// Gaussian fit score: measures how well the XIC intensities across the RT window
        /// match a Gaussian peak shape.
        /// 
        /// Algorithm:
        ///   1. Estimate Gaussian parameters (μ, σ) from weighted moments of the XIC.
        ///      μ = weighted mean RT (weight = intensity)
        ///      σ = weighted standard deviation of RT
        ///   2. Compute predicted Gaussian intensity at each RT point (scaled to match
        ///      the observed apex intensity).
        ///   3. Score = normalized dot product between observed and predicted vectors
        ///      (equivalent to cosine similarity in intensity space).
        /// 
        /// Range [0, 1]; 1 = perfect Gaussian, 0 = orthogonal to any Gaussian.
        /// Returns NaN if fewer than MinMs1Points points or if σ = 0 (all at same RT).
        /// </summary>
        private static float ComputeGaussianFitScore(
            ReadOnlySpan<float> rts,
            ReadOnlySpan<float> intensities)
        {
            if (rts.Length < MinMs1Points) return float.NaN;

            // Step 1: Weighted moments to estimate μ and σ
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

            // Step 2: Build predicted Gaussian vector scaled to observed apex
            // Use stackalloc for small arrays (≤ 256 scans), heap for larger
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

            // Step 3: Cosine similarity between observed and predicted
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

        /// <summary>
        /// Returns the median of all strictly-positive values in the span.
        /// Returns 0f if there are no positive values.
        ///
        /// Retained for reference. Prompt 9 uses <see cref="LowerQuartileNonZero"/>
        /// for PrecursorXicApexIntensity because the median is inflated by the elution
        /// peak body itself when the RT window is wide relative to peptide FWHM.
        ///
        /// For small XIC arrays (≤ 256 scans) the sort buffer is stack-allocated.
        /// Larger arrays fall back to ArrayPool to avoid stack overflow.
        /// </summary>
        private static float MedianNonZero(ReadOnlySpan<float> values)
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

                int mid = count / 2;
                return (count % 2 == 1) ? buf[mid] : (buf[mid - 1] + buf[mid]) * 0.5f;
            }
            finally
            {
                if (rentedBuf != null)
                    ArrayPool<float>.Shared.Return(rentedBuf);
            }
        }

        /// <summary>
        /// Returns the 25th percentile (lower quartile) of all strictly-positive values
        /// in the span. Returns 0f if there are no positive values.
        ///
        /// Used by PrecursorXicApexIntensity (Prompt 9) as the background reference for
        /// log2(apex / Q25) normalization. Q25 represents the genuine low-signal floor
        /// of the XIC window and is not inflated by the elution peak body, unlike the
        /// median. This restores the correct target &gt; decoy direction for this feature.
        ///
        /// Percentile convention: nearest-rank method (index = floor(0.25 * count)).
        /// For small arrays (≤ 256 scans) the sort buffer is stack-allocated.
        /// Larger arrays fall back to ArrayPool to avoid stack overflow.
        /// </summary>
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

                // Nearest-rank: Q25 index = floor(0.25 * count), clamped to valid range.
                // For count=1, this returns buf[0] (the only value, which is the apex itself;
                // log2(apex/Q25) = 0 in that case — no discrimination, which is correct since
                // a single-scan XIC provides no shape information).
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