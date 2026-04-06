// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Represents a single DIA precursor identification result.
    /// This is the DIA analog of a DDA spectral match (PSM), but differs fundamentally:
    /// - Evidence comes from extracted ion chromatograms (XICs) across many scans, not a single spectrum.
    /// - Scores reflect coelution patterns and library-vs-extracted intensity agreement.
    /// 
    /// Lives in MassSpectrometry.Dia (not Omics) because it is tightly coupled to the
    /// DIA engine's SoA extraction output format.
    /// 
    /// Field history:
    /// - Phase 7:  Base fields (DotProductScore, SpectralAngleScore, ExtractedIntensities, XicPointCounts)
    /// - Phase 9-12: Temporal scoring fields (ApexScore, TemporalScore, correlations, peak shape, RT deviation, etc.)
    /// - Phase 13 Prompt 2: Mass accuracy fields (MeanMassErrorPpm, MedianMassErrorPpm, MassErrorStdPpm, MaxAbsMassErrorPpm, ApexObservedMzs)
    /// - Phase 13 Prompt 4: Best-fragment reference curve fields (BestFragCorrelationSum, MedianFragRefCorr, MinFragRefCorr, StdFragRefCorr)
    /// - Phase 13 Prompt 5: Signal ratio deviation fields (MeanSignalRatioDev, MaxSignalRatioDev, StdSignalRatioDev)
    /// - Phase 13 Action 4: Backward-compatible aliases removed (ApexDotProductScore, TemporalCosineScore, etc.)
    /// - Phase 16A Prompt 1: BestFragWeightedCosine, BoundarySignalRatio, ApexToMeanRatio migrated to classifier
    /// - Phase 16B Prompt 5: MS1 feature fields added (PrecursorXicApexIntensity, IsotopePatternScore, Ms1Ms2Correlation, PrecursorElutionScore)
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

        /// <summary>Index into the combined precursors list. Set during assembly. -1 if unknown.</summary>
        public int PrecursorIndex { get; set; } = -1;

        #endregion

        #region Scores

        /// <summary>
        /// Normalized dot product between library and extracted fragment intensities.
        /// Range [0, 1], higher is better. NaN if insufficient fragments.
        /// </summary>
        public float DotProductScore { get; set; }

        /// <summary>
        /// Spectral angle score: 1 - (2/pi) * arccos(normalized dot product).
        /// Range [0, 1], higher is better. NaN if insufficient fragments.
        /// </summary>
        public float SpectralAngleScore { get; set; }

        #endregion

        #region Fragment Evidence

        /// <summary>Number of library fragment ions that yielded at least one XIC data point</summary>
        public int FragmentsDetected { get; set; }

        /// <summary>Total number of library fragment ions queried</summary>
        public int FragmentsQueried { get; }

        /// <summary>
        /// Sum of extracted intensities per fragment (length = FragmentsQueried).
        /// Parallel to the library's fragment ion order.
        /// Used by scorers: library intensities vs these extracted intensities.
        /// </summary>
        public float[] ExtractedIntensities { get; }

        /// <summary>
        /// Number of XIC data points per fragment (length = FragmentsQueried).
        /// A fragment with 0 points was not detected.
        /// </summary>
        public int[] XicPointCounts { get; }

        #endregion

        #region Retention Time Context

        /// <summary>Library/predicted retention time (minutes). Null if library had no RT.</summary>
        public double? LibraryRetentionTime { get; }

        /// <summary>Lower bound of the RT extraction window (minutes)</summary>
        public float RtWindowStart { get; }

        /// <summary>Upper bound of the RT extraction window (minutes)</summary>
        public float RtWindowEnd { get; }

        #endregion



        #region Temporal Scoring (Phase 9-12)

        /// <summary>
        /// Cosine similarity at the apex scan (scan with max summed fragment intensity).
        /// Range [0, 1], higher is better. NaN if apex scan has insufficient signal.
        /// </summary>
        public float ApexScore { get; set; }

        /// <summary>
        /// Temporal cosine score computed from the time x fragment intensity matrix.
        /// Range [0, 1], higher is better.
        /// </summary>
        public float TemporalScore { get; set; }

        /// <summary>
        /// Spectral angle computed from temporal scoring (may differ from SpectralAngleScore
        /// which uses summed intensities).
        /// </summary>
        public float SpectralAngle { get; set; }

        /// <summary>
        /// Mean pairwise Pearson correlation across all detected fragment XICs (full window).
        /// Range [-1, 1], higher is better. NaN if fewer than 3 detected fragments.
        /// </summary>
        public float MeanFragCorr { get; set; }

        /// <summary>
        /// Minimum pairwise Pearson correlation across all detected fragment XICs (full window).
        /// Range [-1, 1]. NaN if fewer than 3 detected fragments.
        /// </summary>
        public float MinFragCorr { get; set; }

        /// <summary>
        /// Mean pairwise Pearson correlation restricted to detected peak group boundaries.
        /// Better than full-window MeanFragCorr because it excludes noise outside the elution peak.
        /// Range [-1, 1], higher is better.
        /// </summary>
        public float PeakMeanFragCorr { get; set; }

        /// <summary>
        /// Minimum pairwise correlation restricted to detected peak boundaries.
        /// Range [-1, 1].
        /// </summary>
        public float PeakMinFragCorr { get; set; }

        /// <summary>
        /// Apex cosine score restricted to peak boundaries.
        /// </summary>
        public float PeakApexScore { get; set; }

        /// <summary>
        /// Peak width in minutes from the detected peak group.
        /// NaN if no peak detected.
        /// </summary>
        public float PeakWidth { get; set; }

        /// <summary>
        /// Peak symmetry ratio (left width / right width). 1.0 = perfectly symmetric.
        /// NaN if no peak detected.
        /// </summary>
        public float PeakSymmetry { get; set; }

        /// <summary>
        /// Number of candidate peak groups detected. Higher values indicate
        /// more ambiguous chromatographic evidence.
        /// </summary>
        public float CandidateCount { get; set; }

        /// <summary>
        /// Log10 of total extracted intensity across all fragments.
        /// </summary>
        public float LogTotalIntensity { get; set; }

        /// <summary>
        /// Coefficient of variation of extracted intensities across detected fragments.
        /// </summary>
        public float IntensityCV { get; set; }

        /// <summary>
        /// Fraction of queried fragments that were detected (had >= 1 XIC data point).
        /// </summary>
        public float FragDetRate { get; set; }

        /// <summary>
        /// Absolute deviation between observed and library retention time (minutes).
        /// </summary>
        public float RtDeviationMinutes { get; set; }

        /// <summary>
        /// Squared deviation between observed and library retention time.
        /// Provides quadratic penalty for larger deviations.
        /// </summary>
        public float RtDeviationSquared { get; set; }

        #endregion

        #region Mass Accuracy (Phase 13, Prompt 2-3)

        /// <summary>
        /// Mean signed mass error in ppm across detected fragments at the apex scan.
        /// Targets cluster near 0; decoys are uniformly distributed across tolerance.
        /// NaN if no fragments matched at the apex scan.
        /// </summary>
        public float MeanMassErrorPpm { get; set; }

        /// <summary>
        /// Median signed mass error in ppm across detected fragments at the apex scan.
        /// More robust to outliers than mean. NaN if no fragments matched.
        /// </summary>
        public float MedianMassErrorPpm { get; set; }

        /// <summary>
        /// Standard deviation of mass errors in ppm across detected fragments.
        /// Targets have tight clustering (low std); decoys have random spread.
        /// NaN if fewer than 2 fragments matched.
        /// </summary>
        public float MassErrorStdPpm { get; set; }

        /// <summary>
        /// Maximum absolute mass error in ppm among detected fragments.
        /// Captures worst-case accuracy. NaN if no fragments matched.
        /// </summary>
        public float MaxAbsMassErrorPpm { get; set; }

        /// <summary>
        /// Observed m/z values at the apex scan for each fragment (parallel to library order).
        /// NaN for fragments not found. Populated by DiaMassAccuracyHelper.
        /// </summary>
        public float[] ApexObservedMzs { get; set; }

        #endregion

        #region Best-Fragment Reference Curve (Phase 13, Prompt 4)

        /// <summary>
        /// Sum of pairwise correlations for the best (least interfered) fragment.
        /// Higher = one fragment has strong coelution with all others.
        /// NaN if fewer than 3 detected fragments.
        /// </summary>
        public float BestFragCorrelationSum { get; set; }

        /// <summary>
        /// Median correlation of each fragment with the best reference fragment.
        /// Range [-1, 1], higher is better.
        /// </summary>
        public float MedianFragRefCorr { get; set; }

        /// <summary>
        /// Minimum correlation with the best reference fragment.
        /// Low values indicate one or more interfered fragments.
        /// </summary>
        public float MinFragRefCorr { get; set; }

        /// <summary>
        /// Standard deviation of correlations with the best reference fragment.
        /// Low values indicate consistent coelution; high values indicate interference.
        /// </summary>
        public float StdFragRefCorr { get; set; }

        #endregion

        #region Signal Ratio Deviation (Phase 13, Prompt 5)

        /// <summary>
        /// Mean deviation of observed fragment intensity ratios from library ratios.
        /// Computed at apex scan. Lower is better for targets.
        /// NaN if fewer than 2 detected fragments at apex.
        /// </summary>
        public float MeanSignalRatioDev { get; set; }

        /// <summary>
        /// Maximum deviation of observed fragment intensity ratios from library ratios.
        /// Captures worst-case interference. Lower is better.
        /// </summary>
        public float MaxSignalRatioDev { get; set; }

        /// <summary>
        /// Standard deviation of signal ratio deviations across fragments.
        /// Lower is better -- consistent ratios indicate clean signal.
        /// </summary>
        public float StdSignalRatioDev { get; set; }

        #endregion

        #region Smoothed Correlations and S/N (Phase 13, Prompt 6)

        /// <summary>
        /// Mean pairwise Pearson correlation on Savitzky-Golay smoothed XICs.
        /// Noise-robust version of MeanFragCorr. Range [-1, 1], higher is better.
        /// NaN if fewer than 2 detected fragments after smoothing.
        /// </summary>
        public float SmoothedMeanFragCorr { get; set; }

        /// <summary>
        /// Minimum pairwise Pearson correlation on smoothed XICs.
        /// Noise-robust version of MinFragCorr. Range [-1, 1], higher is better.
        /// </summary>
        public float SmoothedMinFragCorr { get; set; }

        /// <summary>
        /// Log2 signal-to-noise ratio. Signal = TIC at apex, noise = median TIC across scans.
        /// Higher is better. Typical targets: 2-10. NaN if apex signal is zero.
        /// </summary>
        public float Log2SignalToNoise { get; set; }

        #endregion

        #region Migrated Features (Phase 13 Action Items 1-2)

        /// <summary>
        /// Cosine score at apex weighted by each fragment's correlation with the best fragment.
        /// Fragments that correlate poorly with the reference are down-weighted.
        /// NaN if insufficient fragments.
        /// </summary>
        public float BestFragWeightedCosine { get; set; }

        /// <summary>
        /// Index of the best (least interfered) fragment in the fragment array.
        /// -1 if not computed.
        /// </summary>
        public int BestFragIndex { get; set; }

        /// <summary>
        /// Ratio of boundary signal to apex signal. Low for true chromatographic peaks
        /// (signal drops at boundaries), high for interference (flat signal).
        /// NaN if apex signal is zero.
        /// </summary>
        public float BoundarySignalRatio { get; set; }

        /// <summary>
        /// Ratio of apex TIC to mean TIC across the peak range.
        /// Higher values indicate a sharper, more prominent peak.
        /// NaN if mean signal is zero.
        /// </summary>
        public float ApexToMeanRatio { get; set; }

        #endregion

        #region MS1 Features (Phase 16B, Prompt 5)

        /// <summary>
        /// Log2 of the apex intensity in the precursor (M0) XIC within the calibrated RT window.
        /// High for true matches (precursor elutes with the fragments); low or zero for decoys
        /// (random m/z, no precursor signal at the queried position).
        /// NaN if MS1 scans are absent or the XIC has no signal.
        /// Feature [29].
        /// </summary>
        public float PrecursorXicApexIntensity { get; set; }

        /// <summary>
        /// Dot product of observed [M0, M+1, M+2] intensities at the XIC apex scan vs the
        /// theoretical averagine isotope distribution at the precursor neutral mass.
        /// Range [0, 1]. True peptides follow the averagine pattern; decoy m/z values
        /// land on random background peaks that do not match the expected pattern.
        /// NaN if MS1 scans are absent or fewer than 2 isotopes had signal.
        /// Feature [30].
        /// </summary>
        public float IsotopePatternScore { get; set; }

        /// <summary>
        /// Pearson correlation between the M0 precursor XIC and the best fragment XIC
        /// over the calibrated RT window. True co-eluting signals produce high positive
        /// correlation; noise or interference produces near-zero or negative values.
        /// NaN if MS1 scans are absent or fewer than 3 co-detected time points.
        /// Feature [31].
        /// </summary>
        public float Ms1Ms2Correlation { get; set; }

        /// <summary>
        /// Gaussian-fit score of the precursor XIC peak shape.
        /// Measures how well the M0 XIC across the RT window fits a Gaussian profile.
        /// Range [0, 1]; 1 = perfect Gaussian, 0 = flat or irregular.
        /// Complements ApexScore (which measures fragment co-elution at the apex scan).
        /// NaN if MS1 scans are absent or the XIC has insufficient points.
        /// Feature [32].
        /// </summary>
        public float PrecursorElutionScore { get; set; }

        #endregion

        #region Interference / Chimeric Features (Phase 19, Priority 2)

        /// <summary>
        /// Fraction of this precursor's matched fragment signal that is uncontested —
        /// i.e. not shared with any co-isolated precursor in the same DIA window.
        /// Range [0, 1]. 1.0 = all fragment m/z values are unique to this precursor;
        /// 0.0 = every fragment m/z is also queried by at least one other precursor.
        /// High values indicate a clean, non-chimeric identification.
        /// NaN if not computed (e.g. only one precursor in the window).
        /// Feature [33].
        /// </summary>
        public float ChimericScore { get; set; }

        #endregion

        #region Peak Selection Quality (Prompt 2)

        /// <summary>
        /// Standard deviation of per-fragment apex positions within the detected peak boundaries.
        /// Low value = all fragments co-elute tightly (characteristic of a true peptide peak).
        /// High value = fragments peak at different times (characteristic of interference).
        /// Populated by DiaPeakGroupDetector.SelectBest() via PeakGroupCandidate.CoElutionStd.
        /// 0f if no valid peak detected or only one fragment had signal.
        /// Feature candidate for Prompt 3 feature vector addition.
        /// </summary>
        public float CoElutionStd { get; set; }

        /// <summary>
        /// Difference between the best and second-best candidate SelectionScores
        /// (SelectionScore - SecondBestScore). Large gap = one peak clearly dominates
        /// the selection (confident ID). Small gap = ambiguous selection, two competing peaks.
        /// 0f if only one candidate existed or no valid peak detected.
        /// Feature candidate for Prompt 3 feature vector addition.
        /// </summary>
        public float CandidateScoreGap { get; set; }

        /// <summary>
        /// MS1 apex confirmation score for the selected peak group.
        /// 
        /// The raw ratio of precursor MS1 intensity at the selected apex RT to the maximum
        /// precursor MS1 intensity in the search window. Range [0, 1].
        ///
        ///   1.0 = precursor has strong MS1 signal at the selected apex → confirms target
        ///   0.0 = precursor has no MS1 signal at the selected apex → possible interference
        ///
        /// Default 1.0f (neutral) when MS1 data is unavailable or the search window is
        /// narrow (≤ 1.0 min), in which case MS1 confirmation is not applied.
        ///
        /// This is the raw ratio before sqrt-softening. The sqrt-softened version is used
        /// multiplicatively in SelectBest() to penalize interference candidates; the raw
        /// ratio is stored here as classifier feature [37] for target/decoy separation.
        ///
        /// Feature [37] — Ms1ApexConfirmationScore.
        /// </summary>
        public float Ms1ApexConfirmationScore { get; set; }

        #endregion

        #region Derived RT and Coverage Features (Phase 19, Priority 5)

        /// <summary>
        /// Intensity-weighted fraction of theoretical library fragment ions that were detected.
        /// Each fragment's contribution is weighted by its library intensity, so high-intensity
        /// library fragments that are missing penalize the score more than low-intensity ones.
        /// Range [0, 1]. Higher is better.
        /// Feature [34].
        /// </summary>
        public float LibraryCoverageFraction { get; set; }

        #endregion

        #region FDR and Scoring (set by DiaFdrEngine)

        /// <summary>
        /// Composite classifier score assigned by the FDR engine (LDA or GBT).
        /// Higher scores indicate more confident identifications.
        /// </summary>
        public float ClassifierScore { get; set; }

        /// <summary>
        /// FDR statistics (q-value, cumulative target/decoy counts) assigned by the FDR engine.
        /// Null until FDR estimation is performed.
        /// </summary>
        public DiaFdrInfo FdrInfo { get; set; }

        #endregion

        #region Temporal Metadata (set during assembly)

        /// <summary>
        /// Number of time points that contributed to the temporal cosine score.
        /// </summary>
        public int TimePointsUsed { get; set; }

        /// <summary>
        /// Observed apex retention time (minutes) from the chromatographic peak.
        /// NaN if not determined.
        /// </summary>
        public float ObservedApexRt { get; set; }

        /// <summary>
        /// Detected chromatographic peak group from DiaPeakGroupDetector.SelectBest().
        /// Null if no valid peak was detected.
        /// Prompt 2: PeakGroup struct now contains selection scoring fields
        /// (ApexCosine, CoElutionStd, Snr, SelectionScore, TotalCandidateCount, SecondBestScore).
        /// </summary>
        public PeakGroup? DetectedPeakGroup { get; set; }

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

            // Base scores
            DotProductScore = float.NaN;
            SpectralAngleScore = float.NaN;

            // Temporal scoring (Phase 9-12)
            ApexScore = float.NaN;
            TemporalScore = float.NaN;
            SpectralAngle = float.NaN;
            MeanFragCorr = float.NaN;
            MinFragCorr = float.NaN;
            PeakMeanFragCorr = float.NaN;
            PeakMinFragCorr = float.NaN;
            PeakApexScore = float.NaN;
            PeakWidth = float.NaN;
            PeakSymmetry = float.NaN;
            CandidateCount = float.NaN;
            LogTotalIntensity = float.NaN;
            IntensityCV = float.NaN;
            FragDetRate = float.NaN;
            RtDeviationMinutes = float.NaN;
            RtDeviationSquared = float.NaN;

            // Mass accuracy (Phase 13, Prompt 2-3)
            MeanMassErrorPpm = float.NaN;
            MedianMassErrorPpm = float.NaN;
            MassErrorStdPpm = float.NaN;
            MaxAbsMassErrorPpm = float.NaN;
            ApexObservedMzs = null;

            // Best-fragment reference curve (Phase 13, Prompt 4)
            BestFragCorrelationSum = float.NaN;
            MedianFragRefCorr = float.NaN;
            MinFragRefCorr = float.NaN;
            StdFragRefCorr = float.NaN;

            // Signal ratio deviation (Phase 13, Prompt 5)
            MeanSignalRatioDev = float.NaN;
            MaxSignalRatioDev = float.NaN;
            StdSignalRatioDev = float.NaN;

            // Smoothed correlations and S/N (Phase 13, Prompt 6)
            SmoothedMeanFragCorr = float.NaN;
            SmoothedMinFragCorr = float.NaN;
            Log2SignalToNoise = float.NaN;

            // Migrated features (Phase 13 Action Items 1-2)
            BestFragWeightedCosine = float.NaN;
            BestFragIndex = -1;
            BoundarySignalRatio = float.NaN;
            ApexToMeanRatio = float.NaN;

            // MS1 features (Phase 16B, Prompt 5)
            PrecursorXicApexIntensity = float.NaN;
            IsotopePatternScore = float.NaN;
            Ms1Ms2Correlation = float.NaN;
            PrecursorElutionScore = float.NaN;

            // Interference / chimeric features (Phase 19, Priority 2)
            ChimericScore = float.NaN;

            // Peak selection quality (Prompt 2)
            CoElutionStd = 0f;
            CandidateScoreGap = 0f;

            // MS1 apex confirmation (MS1 Interference Resolution phase)
            // Default 1.0f = neutral (no penalty) until peak group selection computes the real value.
            Ms1ApexConfirmationScore = 1.0f;

            // Derived RT and coverage features (Phase 19, Priority 5)
            // RtDeviationNormalized dropped: 100% NaN (PeakWidth=0 when no peak group detected).
            LibraryCoverageFraction = float.NaN;

            // FDR and scoring
            ClassifierScore = float.NaN;
            FdrInfo = null;

            // Temporal metadata
            TimePointsUsed = 0;
            ObservedApexRt = float.NaN;
            DetectedPeakGroup = null;
        }

        /// <summary>Whether this result meets the minimum fragment detection threshold.</summary>
        public bool MeetsMinFragments(int minRequired) => FragmentsDetected >= minRequired;

        /// <summary>Fraction of queried fragments that were detected (0-1).</summary>
        public float FragmentDetectionRate =>
            FragmentsQueried > 0 ? (float)FragmentsDetected / FragmentsQueried : 0f;

        public override string ToString()
        {
            return $"{Sequence}/{ChargeState} Window={WindowId} " +
                   $"DotProduct={DotProductScore:F4} SpectralAngle={SpectralAngleScore:F4} " +
                   $"Fragments={FragmentsDetected}/{FragmentsQueried}" +
                   (IsDecoy ? " [DECOY]" : "");
        }
    }
}