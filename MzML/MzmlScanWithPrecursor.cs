using MassSpectrometry;
using MzLibUtil;

namespace IO.MzML
{
    public class MzmlScanWithPrecursor : MsDataScanWithPrecursor<MzmlMzSpectrum>, IMzmlScan
    {

        #region Public Constructors

        public MzmlScanWithPrecursor(int ScanNumber, MzmlMzSpectrum massSpectrum, int MsnOrder, bool isCentroid, Polarity Polarity, double RetentionTime, MzRange MzRange, string ScanFilter, MZAnalyzerType MzAnalyzer, double TotalIonCurrent, double selectedIonMz, int? selectedIonChargeStateGuess, double? selectedIonIntensity, double isolationMZ, double? isolationWidth, DissociationType dissociationType, int oneBasedPrecursorScanNumber, double? selectedIonGuessMonoisotopicMZ, double? injectionTime)
            : base(massSpectrum, ScanNumber, MsnOrder, isCentroid, Polarity, RetentionTime, MzRange, ScanFilter, MzAnalyzer, TotalIonCurrent, selectedIonMz, selectedIonChargeStateGuess, selectedIonIntensity, isolationMZ, isolationWidth, dissociationType, oneBasedPrecursorScanNumber, selectedIonGuessMonoisotopicMZ, injectionTime, null)
        {
        }

        #endregion Public Constructors

    }
}