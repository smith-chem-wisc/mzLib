using MassSpectrometry;
using MzLibUtil;

namespace IO.Thermo
{
    public class ThermoScan : MsDataScan<ThermoSpectrum>, IThermoScan
    {

        #region Public Constructors

        public ThermoScan(int oneBasedScanNumber, ThermoSpectrum massSpectrum, int msnOrder, Polarity polarity, double retentionTime, MzRange scanWindowRange, string scanFilter, MZAnalyzerType mzAnalyzer, double totalIonCurrent, double? injectionTime, double[,] noiseData)
        : base(massSpectrum, oneBasedScanNumber, msnOrder, true, polarity, retentionTime, scanWindowRange, scanFilter, mzAnalyzer, totalIonCurrent, injectionTime, noiseData)

        {
        }

        #endregion Public Constructors

    }
}