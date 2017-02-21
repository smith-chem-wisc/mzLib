using MassSpectrometry;
using MzLibUtil;

namespace IO.MzML
{
    public class MzmlScan : MsDataScan<MzmlMzSpectrum>, IMzmlScan
    {

        #region Public Constructors

        public MzmlScan(int oneBasedScanNumber, MzmlMzSpectrum massSpectrum, int msnOrder, bool isCentroid, Polarity polarity, double retentionTime, MzRange scanWindowRange, string scanFilter, MZAnalyzerType mzAnalyzer, double totalIonCurrent)
        : base(oneBasedScanNumber, msnOrder, isCentroid, polarity, retentionTime, scanWindowRange, scanFilter, mzAnalyzer, totalIonCurrent)

        {
            this.MassSpectrum = massSpectrum;
        }

        #endregion Public Constructors

    }
}