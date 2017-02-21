using MassSpectrometry;
using MzLibUtil;

namespace IO.MzML
{
    public class MzmlScanWithPrecursor : MsDataScanWithPrecursor<MzmlMzSpectrum>, IMzmlScan
    {

        #region Public Constructors

        public MzmlScanWithPrecursor(int ScanNumber, MzmlMzSpectrum massSpectrum, int MsnOrder, bool isCentroid, Polarity Polarity, double RetentionTime, MzRange MzRange, string ScanFilter, MZAnalyzerType MzAnalyzer, double TotalIonCurrent, double? selectedIonGuessMZ, int? selectedIonGuessChargeStateGuess, double? selectedIonGuessIntensity, double isolationMZ, double? isolationWidth, DissociationType dissociationType, int oneBasedPrecursorScanNumber, double? selectedIonGuessMonoisotopicMZ, double? injectionTime)
            : base(ScanNumber, MsnOrder, isCentroid, Polarity, RetentionTime, MzRange, ScanFilter, MzAnalyzer, TotalIonCurrent, selectedIonGuessMZ, selectedIonGuessChargeStateGuess, selectedIonGuessIntensity, isolationMZ, isolationWidth, dissociationType, oneBasedPrecursorScanNumber, selectedIonGuessMonoisotopicMZ, injectionTime)
        {
            this.MassSpectrum = massSpectrum;
        }

        #endregion Public Constructors

    }
}