using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class TimsDataScan : MsDataScan
    {
        public int ScanNumberStart { get; }
        public int ScanNumberEnd { get; }
        public double OneOverK0 { get; }
        public int PrecursorId { get; }
        public long FrameId { get; }


        // Need to incorporate scan range somehow
        public TimsDataScan(MzSpectrum massSpectrum,
            int oneBasedScanNumber,
            int msnOrder,
            bool isCentroid,
            Polarity polarity,
            double retentionTime,
            MzRange scanWindowRange,
            string scanFilter,
            MZAnalyzerType mzAnalyzer,
            double totalIonCurrent,
            double? injectionTime,
            double[,] noiseData,
            string nativeId,
            long frameId,
            int scanNumberStart,
            int scanNumberEnd,
            double medianOneOverK0,
            int precursorId,
            double? selectedIonMz = null,
            int? selectedIonChargeStateGuess = null,
            double? selectedIonIntensity = null,
            double? isolationMZ = null,
            double? isolationWidth = null,
            DissociationType? dissociationType = null,
            int? oneBasedPrecursorScanNumber = null,
            double? selectedIonMonoisotopicGuessMz = null,
            string hcdEnergy = null) : 
            base(massSpectrum, oneBasedScanNumber, msnOrder, isCentroid, polarity,
                retentionTime, scanWindowRange, scanFilter, mzAnalyzer, totalIonCurrent,
                injectionTime, noiseData, nativeId, selectedIonMz, selectedIonChargeStateGuess,
                selectedIonIntensity, isolationMZ, isolationWidth, dissociationType,
                oneBasedPrecursorScanNumber, selectedIonMonoisotopicGuessMz, hcdEnergy)
        {
            FrameId = frameId;
            ScanNumberStart = scanNumberStart;
            ScanNumberEnd = scanNumberEnd;
            OneOverK0 = medianOneOverK0;
            PrecursorId = precursorId;
        }
    }
}
