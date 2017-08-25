using MassSpectrometry;
using MzLibUtil;

namespace IO.MzML
{
    public class MzmlScan : MsDataScan<MzmlMzSpectrum>, IMzmlScan
    {
        #region Public Constructors

        public MzmlScan(int oneBasedScanNumber, MzmlMzSpectrum massSpectrum, int msnOrder, bool isCentroid, Polarity polarity, double retentionTime, MzRange scanWindowRange, string scanFilter, MZAnalyzerType mzAnalyzer, double totalIonCurrent, double? injectionTime)
        : base(massSpectrum, oneBasedScanNumber, msnOrder, isCentroid, polarity, retentionTime, scanWindowRange, scanFilter, mzAnalyzer, totalIonCurrent, injectionTime, null)
        {
        }

        #endregion Public Constructors
    }
}