// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Buffers;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Feature vector for a single DIA precursor identification.
    /// 
    /// Phase 13 Action Item 5: 29-feature classifier vector.
    /// All 3 migrated features now active (BestFragWeightedCosine, BoundarySignalRatio,
    /// remain commented out until Prompt 5 activates them (26 -> 29 features).
    /// 
    /// DROPPED (from Phase 12/early Phase 13):
    ///   - MeanFragCorr / MinFragCorr -- replaced by SmoothedMeanFragCorr/SmoothedMinFragCorr + best-fragment features
    ///   - PeakApexScore -- negative separation in Phase 12 analysis (hurts discrimination)
    ///   - PeakSymmetry -- zero separation (0.01)
    ///   - MedianXicDepth / XicDepthCV -- weak features, not in Prompt 8 plan
    /// 
    /// ADDED:
    ///   - CandidateCount (strong discriminator, separation 1.66)
    ///   - MeanMassErrorPpm (absolute), MassErrorStdPpm, MaxAbsMassErrorPpm (mass accuracy: 3)
    ///   - SmoothedMinFragCorr (was computed but not in vector)
    /// 
    /// Current layout: 29 classifier features.
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

        // -- Metadata (not classifier features) ---------------------------
        public bool IsDecoy;
        public int PrecursorIndex;
        public int ChargeState;
        public float PrecursorMz;
        public int FragmentsDetected;
        public int FragmentsQueried;

        /// <summary>
        /// Number of features used by the classifier.
        /// Phase 13 Action Item 5: 26 -> 29 (activated 3 migrated features).
        /// </summary>
        public const int ClassifierFeatureCount = 29;

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
        /// Phase 13 Action Item 5 order (29 features):
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
        };
    }

    /// <summary>
    /// Computes feature vectors from DiaSearchResult objects.
    /// Thread-safe: uses only local state + ArrayPool.
    /// 
    /// Phase 13 Action Item 5: 29-feature vector with all migrated features active.
    /// </summary>
    public static class DiaFeatureExtractor
    {
        /// <summary>
        /// Computes a feature vector from a scored DiaSearchResult.
        /// </summary>
        public static DiaFeatureVector ComputeFeatures(
            DiaSearchResult result,
            int precursorIndex)
        {
            var fv = new DiaFeatureVector();

            // -- Primary scores [0-2] ------------------------------------
            fv.ApexScore = SafeScore(result.ApexScore);
            fv.TemporalScore = SafeScore(result.TemporalScore);
            fv.SpectralAngle = SafeScore(result.SpectralAngleScore);

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
    }
}