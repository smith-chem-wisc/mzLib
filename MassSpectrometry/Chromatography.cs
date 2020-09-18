using Chemistry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MassSpectrometry
{
    public class Chromatography
    {
        /// <summary>
        /// Extracts an ion chromatogram from a spectra file, given a mass, charge, retention time, and mass tolerance.
        /// </summary>
        public static ExtractedIonChromatogram ExtractIonChromatogram(MsDataFile file, double mass, int charge, Tolerance massTolerance,
            double retentionTime, int msOrder = 1, double retentionTimeWindowWidth = 5)
        {
            double theorMz = mass.ToMz(charge);
            double startRt = retentionTime - retentionTimeWindowWidth / 2;
            double endRt = retentionTime + retentionTimeWindowWidth / 2;
            List<Datum> xicData = new List<Datum>();

            IEnumerable<MsDataScan> scansInRtWindow = file.GetMsScansInTimeRange(startRt, endRt);

            foreach (MsDataScan scan in scansInRtWindow.Where(p => p.MsnOrder == msOrder))
            {
                int ind = scan.MassSpectrum.GetClosestPeakIndex(theorMz);

                double expMz = scan.MassSpectrum.XArray[ind];

                if (massTolerance.Within(expMz.ToMass(charge), mass))
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
