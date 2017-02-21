using MassSpectrometry;
using MzLibUtil;

namespace IO.Thermo
{
    public class ThermoScanWithPrecursor : MsDataScanWithPrecursor<ThermoSpectrum>, IThermoScan
    {

        #region Public Constructors

        public ThermoScanWithPrecursor(int ScanNumber, ThermoSpectrum massSpectrum, int MsnOrder, Polarity Polarity, double RetentionTime, MzRange MzRange, string ScanFilter, MZAnalyzerType MzAnalyzer, double TotalIonCurrent, double? selectedIonGuessMZ, int? selectedIonGuessChargeStateGuess, double isolationMZ, double? isolationWidth, DissociationType dissociationType, int oneBasedPrecursorScanNumber, double? selectedIonGuessMonoisotopicMZ)
            : base(ScanNumber, MsnOrder, true, Polarity, RetentionTime, MzRange, ScanFilter, MzAnalyzer, TotalIonCurrent, selectedIonGuessMZ, selectedIonGuessChargeStateGuess, null, isolationMZ, isolationWidth, dissociationType, oneBasedPrecursorScanNumber, selectedIonGuessMonoisotopicMZ)
        {
            this.MassSpectrum = massSpectrum;
        }

        #endregion Public Constructors

    }
}