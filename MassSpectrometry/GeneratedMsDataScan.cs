using MzLibUtil;

namespace MassSpectrometry
{
    internal class GeneratedMsDataScan : MsDataScan<IMzSpectrum<IMzPeak>>
    {

        #region Public Constructors

        public GeneratedMsDataScan(IMzSpectrum<IMzPeak> massSpectrum, int oneBasedScanNumber, int msnOrder, bool isCentroid, Polarity polarity, double retentionTime, MzRange scanWindowRange, string scanFilter, MZAnalyzerType mzAnalyzer, double totalIonCurrent, double? injectionTime, double[,] noiseData)
            : base(oneBasedScanNumber, msnOrder, isCentroid, polarity, retentionTime, scanWindowRange, scanFilter, mzAnalyzer, totalIonCurrent, injectionTime, noiseData)
        {
            MassSpectrum = massSpectrum;
        }

        #endregion Public Constructors

    }
}