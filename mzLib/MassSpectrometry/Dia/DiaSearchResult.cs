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
        /// Normalized dot product between library and extracted fragment intensities.
        /// Range [0, 1], higher is better. NaN if insufficient fragments.
        /// </summary>
        public float DotProductScore { get; set; }

        /// <summary>
        /// Spectral angle score: 1 - (2/π) * arccos(normalized dot product).
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

        #region FDR

        /// <summary>
        /// FDR statistics for this result (q-value, cumulative counts, PEP).
        /// Populated by FdrAnalysisEngineDia during post-search analysis.
        /// Null until FDR analysis is run.
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
        }

        /// <summary>Whether this result meets the minimum fragment detection threshold.</summary>
        public bool MeetsMinFragments(int minRequired) => FragmentsDetected >= minRequired;

        /// <summary>Fraction of queried fragments that were detected (0–1).</summary>
        public float FragmentDetectionRate =>
            FragmentsQueried > 0 ? (float)FragmentsDetected / FragmentsQueried : 0f;

        public override string ToString()
        {
            string fdrStr = FdrInfo != null ? $" q={FdrInfo.QValue:F4}" : "";
            return $"{Sequence}/{ChargeState} Window={WindowId} " +
                   $"DotProduct={DotProductScore:F4} SpectralAngle={SpectralAngleScore:F4} " +
                   $"Fragments={FragmentsDetected}/{FragmentsQueried}" +
                   fdrStr +
                   (IsDecoy ? " [DECOY]" : "");
        }
    }
}
