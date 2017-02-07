using MassSpectrometry;
using MzLibUtil;

namespace IO.MzML
{
    public class MzmlScan : MsDataScan<MzmlMzSpectrum, MzmlPeak>, IMzmlScan
    {
        public MzmlScan(int oneBasedScanNumber, MzmlMzSpectrum massSpectrum, string id, int msnOrder, bool isCentroid, Polarity polarity, double retentionTime, MzRange scanWindowRange, string scanFilter, MZAnalyzerType mzAnalyzer, double injectionTime, double totalIonCurrent)
        : base(oneBasedScanNumber, id, msnOrder, isCentroid, polarity, retentionTime, scanWindowRange, scanFilter, mzAnalyzer, injectionTime, totalIonCurrent)

        {
            this.MassSpectrum = massSpectrum;
        }
    }
}