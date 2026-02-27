// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// FDR statistics for a single DIA precursor match.
    /// 
    /// Simpler than the DDA FdrInfo (EngineLayer.FdrAnalysis.FdrInfo) because DIA results
    /// have no notch-level ambiguity and no multi-hypothesis matches — each DiaSearchResult
    /// maps to exactly one target or decoy library entry.
    /// 
    /// Lives in MassSpectrometry.Dia alongside DiaSearchResult to avoid circular dependencies
    /// between mzLib and MetaMorpheus EngineLayer.
    /// </summary>
    public class DiaFdrInfo
    {
        /// <summary>Cumulative target count at this result's rank position</summary>
        public double CumulativeTarget { get; set; }

        /// <summary>Cumulative decoy count at this result's rank position</summary>
        public double CumulativeDecoy { get; set; }

        /// <summary>
        /// Global q-value: the minimum FDR at which this result would be accepted.
        /// Calculated as cumulative_decoy / max(cumulative_target, 1), then monotonized
        /// (walking from worst to best score, taking the running minimum).
        /// Range [0, 1]. Default 2.0 indicates "not yet calculated".
        /// </summary>
        public double QValue { get; set; }

        /// <summary>
        /// Peptide-level q-value: FDR calculated after collapsing to unique sequences
        /// (best score per sequence). Null if peptide-level FDR has not been run.
        /// </summary>
        public double? PeptideQValue { get; set; }

        /// <summary>
        /// Posterior error probability from a future discriminative model.
        /// NaN until PEP analysis is implemented.
        /// </summary>
        public double PEP { get; set; }

        /// <summary>
        /// Q-value computed from PEP-sorted results.
        /// NaN until PEP analysis is implemented.
        /// </summary>
        public double PEP_QValue { get; set; }

        public DiaFdrInfo()
        {
            QValue = 2.0;       // sentinel: not yet calculated (matches DDA FdrInfo convention)
            PeptideQValue = null;
            PEP = double.NaN;
            PEP_QValue = double.NaN;
        }
    }
}
