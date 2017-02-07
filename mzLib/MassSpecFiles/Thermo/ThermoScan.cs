using MassSpectrometry;
using MzLibUtil;

namespace IO.Thermo
{
    public class ThermoScan : MsDataScan<ThermoSpectrum, ThermoMzPeak>, IThermoScan
    {
        #region Public Constructors

        public ThermoScan(int oneBasedScanNumber, ThermoSpectrum massSpectrum, string id, int msnOrder, bool isCentroid, Polarity polarity, double retentionTime, MzRange scanWindowRange, string scanFilter, MZAnalyzerType mzAnalyzer, double injectionTime, double totalIonCurrent)
        : base(oneBasedScanNumber, id, msnOrder, isCentroid, polarity, retentionTime, scanWindowRange, scanFilter, mzAnalyzer, injectionTime, totalIonCurrent)

        {
            this.MassSpectrum = massSpectrum;
        }

        #endregion Public Constructors
    }
}