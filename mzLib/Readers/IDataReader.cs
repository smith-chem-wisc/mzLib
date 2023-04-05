using Chemistry;
using MzLibUtil;
using MassSpectrometry; 
namespace Readers
{
    public static class MsDataFileHelpers
    {
        
        /// <summary>
        /// Extracts an ion chromatogram from the spectra file, given a mass, charge, retention time, and mass tolerance.
        /// </summary>
        public static ExtractedIonChromatogram ExtractIonChromatogram(this MsDataFile file, double neutralMass, int charge, Tolerance massTolerance,
            double retentionTimeInMinutes, int msOrder = 1, double retentionTimeWindowWidthInMinutes = 5)
        {
            double theorMz = neutralMass.ToMz(charge);
            double startRt = retentionTimeInMinutes - retentionTimeWindowWidthInMinutes / 2;
            double endRt = retentionTimeInMinutes + retentionTimeWindowWidthInMinutes / 2;
            List<Datum> xicData = new List<Datum>();

            IEnumerable<MsDataScan> scansInRtWindow = file.GetMsScansInTimeRange(startRt, endRt);

            foreach (MsDataScan scan in scansInRtWindow.Where(p => p.MsnOrder == msOrder))
            {
                int ind = scan.MassSpectrum.GetClosestPeakIndex(theorMz);

                double expMz = scan.MassSpectrum.XArray[ind];

                if (massTolerance.Within(expMz.ToMass(charge), neutralMass))
                {
                    xicData.Add(new Datum(scan.RetentionTime, scan.MassSpectrum.YArray[ind]));
                }
                else
                {
                    xicData.Add(new Datum(scan.RetentionTime, 0));
                }
            }

            return new ExtractedIonChromatogram(xicData);
        }

    }
}