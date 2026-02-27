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

        #region Scores

        /// <summary>
        /// Primary similarity score between library and extracted fragment intensities.
        /// Range [0, 1], higher is better. NaN if insufficient fragments.
        /// 
        /// The interpretation depends on the ScoringStrategy used:
        ///   Summed: L2-normalized dot product of summed intensity vectors.
        ///   ConsensusApex: cosine at the chromatographic peak.
        ///   TemporalCosine: average cosine across RT time points.
        ///   WeightedTemporalCosineWithTransform: (weighted average cosine)^N.
        /// </summary>
        public float DotProductScore { get; set; }

        /// <summary>
        /// Spectral angle score: 1 - (2/pi) * arccos(normalized dot product).
        /// Range [0, 1], higher is better. NaN if insufficient fragments.
        /// Computed from DotProductScore (or RawCosine when available).
        /// </summary>
        public float SpectralAngleScore { get; set; }

        /// <summary>
        /// The cosine similarity before any nonlinear transform was applied.
        /// Equals DotProductScore for Summed, ConsensusApex, and TemporalCosine strategies.
        /// For WeightedTemporalCosineWithTransform, this is the pre-transform score.
        /// Useful for diagnostics and multi-feature classifiers.
        /// </summary>
        public float RawCosine { get; set; }

        #endregion

        #region Hybrid Scoring Features

        /// <summary>
        /// Cosine similarity at the consensus chromatographic apex (single best time point).
        /// Always computed regardless of which ScoringStrategy was used for DotProductScore.
        /// 
        /// High apex score = the library pattern is well-reproduced at peak intensity.
        /// Good discriminator for high-quality matches (sharp peaks), but noisy for 
        /// low-abundance or interfered peptides.
        /// Range [0, 1], NaN if insufficient data.
        /// </summary>
        public float ApexDotProductScore { get; set; }

        /// <summary>
        /// Average cosine similarity across all RT time points with sufficient fragment signal.
        /// Always computed regardless of which ScoringStrategy was used for DotProductScore.
        /// 
        /// High temporal score = the library pattern is consistently reproduced across the peak.
        /// More robust than apex (averages out noise), but can be pulled down by co-elution.
        /// Range [0, 1], NaN if insufficient data.
        /// </summary>
        public float TemporalCosineScore { get; set; }

        #endregion

        #region Temporal Scoring Diagnostics

        /// <summary>
        /// Number of RT time points that contributed to the temporal score.
        /// Higher values mean more temporal evidence and more robust scoring.
        /// 0 for Summed scoring (not time-resolved).
        /// </summary>
        public int TimePointsUsed { get; set; }

        /// <summary>
        /// Index of the consensus apex time point within the RT window.
        /// The time point where total fragment signal was highest.
        /// -1 if not applicable (e.g., Summed scoring).
        /// </summary>
        public int ApexTimeIndex { get; set; }

        /// <summary>
        /// The scoring strategy that was used to compute DotProductScore.
        /// Stored here so downstream code knows how to interpret the scores.
        /// </summary>
        public ScoringStrategy ScoringStrategyUsed { get; set; }

        #endregion

        #region Fragment Evidence

        /// <summary>Number of library fragment ions that yielded at least one XIC data point</summary>
        public int FragmentsDetected { get; set; }

        /// <summary>Total number of library fragment ions queried</summary>
        public int FragmentsQueried { get; }

        /// <summary>
        /// Sum of extracted intensities per fragment (length = FragmentsQueried).
        /// Parallel to the library's fragment ion order.
        /// Retained for backward compatibility and diagnostic use.
        /// Note: when temporal scoring is used, DotProductScore is NOT computed from
        /// these summed values -- it uses RT-resolved data instead.
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

        /// <summary>
        /// Observed retention time at the chromatographic apex (minutes).
        /// Computed during temporal scoring from the RT of the scan with
        /// maximum total fragment signal.
        /// NaN if no apex could be determined.
        /// </summary>
        public float ObservedApexRt { get; set; }

        #endregion

        #region Coelution Features

        /// <summary>
        /// Mean pairwise Pearson correlation across all detected fragment XICs.
        /// Range [-1, 1], higher is better. True coeluting fragments correlate ~0.9+.
        /// NaN if fewer than 2 fragments have sufficient XIC data points.
        /// 
        /// This is among the most important features for distinguishing true
        /// identifications from interference, per DIA-NN's feature engineering.
        /// </summary>
        public float MeanFragmentCorrelation { get; set; }

        /// <summary>
        /// Minimum pairwise Pearson correlation across all detected fragment XICs.
        /// Identifies the "weakest link" -- the fragment most likely to be interfered.
        /// NaN if fewer than 2 fragments have sufficient XIC data points.
        /// </summary>
        public float MinFragmentCorrelation { get; set; }

        #endregion
        #region Classifier & FDR (Phase 11)

        /// <summary>
        /// Composite classifier score from DiaLinearDiscriminant.
        /// Set by the FDR engine during iterative retraining.
        /// Higher is better. Default NaN = not yet scored.
        /// </summary>
        public float ClassifierScore { get; set; }

        /// <summary>
        /// FDR statistics: q-value, cumulative targets/decoys.
        /// Null until DiaFdrEngine.ComputeQValues() is called.
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
            TimePointsUsed = 0;
            ApexTimeIndex = -1;
            ScoringStrategyUsed = ScoringStrategy.Summed;
            ObservedApexRt = float.NaN;
            MeanFragmentCorrelation = float.NaN;
            MinFragmentCorrelation = float.NaN;
        }

        /// <summary>Whether this result meets the minimum fragment detection threshold.</summary>
        public bool MeetsMinFragments(int minRequired) => FragmentsDetected >= minRequired;

        /// <summary>Fraction of queried fragments that were detected (0-1).</summary>
        public float FragmentDetectionRate =>
            FragmentsQueried > 0 ? (float)FragmentsDetected / FragmentsQueried : 0f;

        public override string ToString()
        {
            string strategyLabel = ScoringStrategyUsed switch
            {
                ScoringStrategy.Summed => "sum",
                ScoringStrategy.ConsensusApex => "apex",
                ScoringStrategy.TemporalCosine => "temporal",
                ScoringStrategy.WeightedTemporalCosineWithTransform => "wt-temporal",
                _ => "?"
            };
            return $"{Sequence}/{ChargeState} Window={WindowId} " +
                   $"DP={DotProductScore:F4}({strategyLabel}) Apex={ApexDotProductScore:F4} Temporal={TemporalCosineScore:F4} " +
                   $"Corr={MeanFragmentCorrelation:F3} " +
                   $"Fragments={FragmentsDetected}/{FragmentsQueried} " +
                   $"TimePts={TimePointsUsed}" +
                   (IsDecoy ? " [DECOY]" : "");
        }
    }
}