// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Represents a single DIA precursor identification result with multi-feature scoring.
    /// 
    /// Phase 13 expansion: from 17 features to ~28 features for improved target/decoy
    /// discrimination. Features are organized into categories matching DIA-NN's architecture.
    /// 
    /// NOTE: ScoringStrategy enum lives in DiaTemporalScorer.cs — do NOT duplicate here.
    /// </summary>
    public class DiaSearchResult
    {
        #region Identification

        /// <summary>Peptide sequence with modifications, as stored in LibrarySpectrum</summary>
        public string Sequence { get; }

        /// <summary>Precursor charge state</summary>
        public int ChargeState { get; }

        /// <summary>Precursor m/z (from library)</summary>
        public double PrecursorMz { get; }

        /// <summary>DIA isolation window ID this precursor was extracted from</summary>
        public int WindowId { get; }

        /// <summary>Whether this is a decoy identification (from decoy library)</summary>
        public bool IsDecoy { get; }

        #endregion

        #region Category 1: Spectral Similarity Scores

        /// <summary>
        /// Normalized dot product between library and extracted fragment intensities.
        /// Range [0, 1], higher is better. NaN if insufficient fragments.
        /// </summary>
        public float DotProductScore { get; set; }

        /// <summary>
        /// Spectral angle score: 1 - (2/pi) * arccos(normalized dot product).
        /// Range [0, 1], higher is better.
        /// </summary>
        public float SpectralAngleScore { get; set; }

        /// <summary>Raw cosine value before clamping/transform.</summary>
        public float RawCosine { get; set; }

        /// <summary>Cosine similarity at the apex time point (full-window apex).</summary>
        public float ApexDotProductScore { get; set; }

        /// <summary>Intensity-weighted average cosine across all time points in the window.</summary>
        public float TemporalCosineScore { get; set; }

        #endregion

        #region Category 2: Fragment Co-elution / Correlation

        /// <summary>Mean pairwise Pearson correlation between fragment XICs (full window).</summary>
        public float MeanFragmentCorrelation { get; set; }

        /// <summary>Minimum pairwise Pearson correlation between fragment XICs (full window).</summary>
        public float MinFragmentCorrelation { get; set; }

        /// <summary>Mean pairwise Pearson correlation within detected peak boundaries.</summary>
        public float PeakMeanFragCorrelation { get; set; }

        /// <summary>Minimum pairwise Pearson correlation within detected peak boundaries.</summary>
        public float PeakMinFragCorrelation { get; set; }

        // ── Phase 13: Best-Fragment Reference Curve ──

        /// <summary>Sum of Pearson correlations between the best fragment and all others.</summary>
        public float BestFragCorrelationSum { get; set; }

        /// <summary>Median correlation of individual fragments with the best reference fragment.</summary>
        public float MedianFragRefCorr { get; set; }

        /// <summary>Minimum correlation of any fragment with the best reference fragment.</summary>
        public float MinFragRefCorr { get; set; }

        /// <summary>Std deviation of fragment-vs-reference correlations.</summary>
        public float StdFragRefCorr { get; set; }

        /// <summary>Cosine at apex weighted by each fragment's correlation with the best reference.</summary>
        public float BestFragWeightedCosine { get; set; }

        /// <summary>Index of the best (least-interfered) fragment.</summary>
        public int BestFragIndex { get; set; }

        // ── Phase 13: Smoothed Correlations ──

        /// <summary>Mean pairwise Pearson correlation on min-of-3-consecutive smoothed XICs.</summary>
        public float SmoothedMeanFragCorr { get; set; }

        /// <summary>Minimum pairwise correlation on smoothed XICs.</summary>
        public float SmoothedMinFragCorr { get; set; }

        #endregion

        #region Category 3: Signal/Noise and Intensity

        /// <summary>Log2 of total extracted intensity within the peak/window.</summary>
        public float LogTotalIntensity { get; set; }

        /// <summary>Log2(signal / noise).</summary>
        public float Log2SignalToNoise { get; set; }

        // ── Phase 13: Per-Fragment Signal Ratio Deviation ──

        /// <summary>Mean |log2(obs_fraction / lib_fraction)| per fragment. LOWER = better.</summary>
        public float MeanSignalRatioDeviation { get; set; }

        /// <summary>Max |log2(obs_fraction / lib_fraction)| across fragments.</summary>
        public float MaxSignalRatioDeviation { get; set; }

        /// <summary>Std of per-fragment signal ratio deviations.</summary>
        public float StdSignalRatioDeviation { get; set; }

        #endregion

        #region Category 4: Peak Shape

        /// <summary>Cosine at the peak-detected apex.</summary>
        public float PeakApexScore { get; set; }

        /// <summary>Temporal cosine restricted to detected peak boundaries.</summary>
        public float PeakTemporalScore { get; set; }

        /// <summary>Ratio of average boundary signal to apex signal. LOWER = better.</summary>
        public float BoundarySignalRatio { get; set; }

        /// <summary>Ratio of apex signal to mean signal. HIGHER = better.</summary>
        public float ApexToMeanRatio { get; set; }

        /// <summary>Number of candidate peak groups detected.</summary>
        public int CandidateCount { get; set; }

        /// <summary>Detected peak group from DiaPeakGroupDetector. Null if none detected.</summary>
        public PeakGroup? DetectedPeakGroup { get; set; }

        /// <summary>Peak width in minutes.</summary>
        public float PeakWidth { get; set; }

        #endregion

        #region Category 5: Retention Time Quality

        /// <summary>Observed apex RT from extraction (minutes).</summary>
        public float ObservedApexRt { get; set; }

        /// <summary>Time point index of the observed apex.</summary>
        public int ApexTimeIndex { get; set; }

        /// <summary>Scoring strategy actually used for this result.</summary>
        public ScoringStrategy ScoringStrategyUsed { get; set; }

        #endregion

        #region Fragment Evidence

        /// <summary>Number of library fragment ions that yielded at least one XIC data point</summary>
        public int FragmentsDetected { get; set; }

        /// <summary>Total number of library fragment ions queried</summary>
        public int FragmentsQueried { get; }

        /// <summary>Sum of extracted intensities per fragment (length = FragmentsQueried).</summary>
        public float[] ExtractedIntensities { get; }

        /// <summary>Number of XIC data points per fragment (length = FragmentsQueried).</summary>
        public int[] XicPointCounts { get; }

        /// <summary>Number of time points used for temporal scoring.</summary>
        public int TimePointsUsed { get; set; }

        #endregion

        #region Retention Time Context

        /// <summary>Library/predicted retention time (minutes). Null if library had no RT.</summary>
        public double? LibraryRetentionTime { get; }

        /// <summary>Lower bound of the RT extraction window (minutes)</summary>
        public float RtWindowStart { get; }

        /// <summary>Upper bound of the RT extraction window (minutes)</summary>
        public float RtWindowEnd { get; }

        #endregion

        #region Classifier and FDR

        /// <summary>
        /// Combined score from the linear discriminant / logistic regression classifier.
        /// Set by DiaLinearDiscriminant after feature extraction. Used for ranking and FDR.
        /// </summary>
        public float ClassifierScore { get; set; }

        /// <summary>
        /// FDR information computed by DiaFdrEngine (q-value, cumulative counts, etc.).
        /// </summary>
        public DiaFdrInfo FdrInfo { get; set; }

        #endregion

        public DiaSearchResult(
            string sequence,
            int chargeState,
            double precursorMz,
            int windowId,
            bool isDecoy,
            int fragmentsQueried,
            double? libraryRetentionTime,
            float rtWindowStart,
            float rtWindowEnd)
        {
            Sequence = sequence ?? throw new ArgumentNullException(nameof(sequence));
            ChargeState = chargeState;
            PrecursorMz = precursorMz;
            WindowId = windowId;
            IsDecoy = isDecoy;
            FragmentsQueried = fragmentsQueried;
            ExtractedIntensities = new float[fragmentsQueried];
            XicPointCounts = new int[fragmentsQueried];
            LibraryRetentionTime = libraryRetentionTime;
            RtWindowStart = rtWindowStart;
            RtWindowEnd = rtWindowEnd;

            DotProductScore = float.NaN;
            SpectralAngleScore = float.NaN;
            RawCosine = float.NaN;
            ApexDotProductScore = float.NaN;
            TemporalCosineScore = float.NaN;
            MeanFragmentCorrelation = float.NaN;
            MinFragmentCorrelation = float.NaN;
            PeakMeanFragCorrelation = float.NaN;
            PeakMinFragCorrelation = float.NaN;
            BestFragCorrelationSum = float.NaN;
            MedianFragRefCorr = float.NaN;
            MinFragRefCorr = float.NaN;
            StdFragRefCorr = float.NaN;
            BestFragWeightedCosine = float.NaN;
            BestFragIndex = -1;
            SmoothedMeanFragCorr = float.NaN;
            SmoothedMinFragCorr = float.NaN;
            LogTotalIntensity = 0f;
            Log2SignalToNoise = float.NaN;
            MeanSignalRatioDeviation = float.NaN;
            MaxSignalRatioDeviation = float.NaN;
            StdSignalRatioDeviation = float.NaN;
            PeakApexScore = float.NaN;
            PeakTemporalScore = float.NaN;
            BoundarySignalRatio = float.NaN;
            ApexToMeanRatio = float.NaN;
            CandidateCount = 0;
            PeakWidth = 0f;
            ClassifierScore = float.NaN;
            FdrInfo = null;
        }

        /// <summary>Whether this result meets the minimum fragment detection threshold.</summary>
        public bool MeetsMinFragments(int minRequired) => FragmentsDetected >= minRequired;

        /// <summary>Fraction of queried fragments that were detected (0-1).</summary>
        public float FragmentDetectionRate =>
            FragmentsQueried > 0 ? (float)FragmentsDetected / FragmentsQueried : 0f;

        /// <summary>RT deviation from predicted in minutes (signed).</summary>
        public float RtDeviationMinutes =>
            LibraryRetentionTime.HasValue ? ObservedApexRt - (float)LibraryRetentionTime.Value : 0f;

        /// <summary>Squared RT deviation (always positive).</summary>
        public float RtDeviationSquared => RtDeviationMinutes * RtDeviationMinutes;

        /// <summary>Intensity coefficient of variation across fragments at apex.</summary>
        public float IntensityCV
        {
            get
            {
                if (FragmentsDetected < 2) return float.NaN;
                float sum = 0f, sum2 = 0f;
                int n = 0;
                for (int i = 0; i < FragmentsQueried; i++)
                {
                    if (ExtractedIntensities[i] > 0f)
                    {
                        sum += ExtractedIntensities[i];
                        sum2 += ExtractedIntensities[i] * ExtractedIntensities[i];
                        n++;
                    }
                }
                if (n < 2) return float.NaN;
                float mean = sum / n;
                float variance = (sum2 / n) - (mean * mean);
                return mean > 0f ? MathF.Sqrt(Math.Max(0f, variance)) / mean : float.NaN;
            }
        }

        public override string ToString()
        {
            return $"{Sequence}/{ChargeState} Window={WindowId} " +
                   $"DotProduct={DotProductScore:F4} SpectralAngle={SpectralAngleScore:F4} " +
                   $"Fragments={FragmentsDetected}/{FragmentsQueried}" +
                   (IsDecoy ? " [DECOY]" : "");
        }
    }
}