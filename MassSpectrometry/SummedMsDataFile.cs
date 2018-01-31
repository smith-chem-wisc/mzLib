using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    public class SummedMsDataFile : MsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>>
    {
        #region Private Fields

        private readonly IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> raw;
        private readonly int numScansToAverage;
        private readonly double ppmToleranceForPeakCombination;

        #endregion Private Fields

        #region Public Constructors

        public SummedMsDataFile(IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> raw, int numScansToAverage, double ppmToleranceForPeakCombination) :
            base(raw.NumSpectra - numScansToAverage + 1,
                new SourceFile(
                @"scan number only nativeID format",
                raw.SourceFile.MassSpectrometerFileFormat,
                raw.SourceFile.CheckSum,
                raw.SourceFile.FileChecksumType,
                raw.SourceFile.Uri,
                raw.SourceFile.Id,
                raw.SourceFile.FileName))
        {
            this.raw = raw;
            this.numScansToAverage = numScansToAverage;
            this.ppmToleranceForPeakCombination = ppmToleranceForPeakCombination;
        }

        #endregion Public Constructors

        #region Public Methods

        public override IEnumerable<IMsDataScan<IMzSpectrum<IMzPeak>>> GetMS1Scans()
        {
            for (int i = 1; i <= NumSpectra; i++)
            {
                var scan = GetOneBasedScan(i);
                if (scan.MsnOrder == 1)
                    yield return scan;
            }
        }

        public override IMsDataScan<IMzSpectrum<IMzPeak>> GetOneBasedScan(int oneBasedScanNumber)
        {
            if (Scans[oneBasedScanNumber - 1] == null)
            {
                var representativeScanNumber = oneBasedScanNumber + (numScansToAverage - 1) / 2;
                var representative = raw.GetOneBasedScan(representativeScanNumber);
                if (representative.MsnOrder != 1)
                    throw new MzLibException("Scan " + representativeScanNumber + " is not MS1 scan");
                int msnOrder = 1;
                Polarity polarity = representative.Polarity;
                if (!representative.IsCentroid)
                    throw new MzLibException("Scan " + representativeScanNumber + " is not centroid scan");
                bool isCentroid = true;
                double retentionTime = representative.RetentionTime;
                MZAnalyzerType mzAnalyzer = representative.MzAnalyzer;

                IMzSpectrum<IMzPeak> peaks = CombinePeaks(raw.Where(b => b.OneBasedScanNumber >= oneBasedScanNumber && b.OneBasedScanNumber <= oneBasedScanNumber + numScansToAverage - 1).Select(b => b.MassSpectrum).ToList(), ppmToleranceForPeakCombination);

                MzRange scanWindowRange = representative.ScanWindowRange;

                double totalIonCurrent = peaks.SumOfAllY;
                double injectionTime = double.NaN;
                double[,] noiseData = null;

                Scans[oneBasedScanNumber - 1] = new MsDataScan<IMzSpectrum<IMzPeak>>(peaks, oneBasedScanNumber, msnOrder, isCentroid, polarity, retentionTime, scanWindowRange, null, mzAnalyzer, totalIonCurrent, injectionTime, noiseData, "scan=" + oneBasedScanNumber);
            }
            return Scans[oneBasedScanNumber - 1];
        }

        #endregion Public Methods

        #region Private Methods

        private IMzSpectrum<IMzPeak> CombinePeaks(List<IMzSpectrum<IMzPeak>> spectraToCombine, double ppmTolerance)
        {
            List<IMzPeak> finalizedPeaks = new List<IMzPeak>();

            int[] peaksLeft = spectraToCombine.Select(b => b.Size).ToArray();
            int[] totalPeaks = spectraToCombine.Select(b => b.Size).ToArray();
            double[] nextPeakMzs = spectraToCombine.Select(b => b.XArray[0]).ToArray();
            double[] nextPeaksIntensites = spectraToCombine.Select(b => b.YArray[0]).ToArray();

            double nextMz = nextPeakMzs.Min();
            int indexOfNextScanToConsider = Array.IndexOf(nextPeakMzs, nextMz);
            GeneratedPeak lastPeak = new GeneratedPeak(nextMz, nextPeaksIntensites[indexOfNextScanToConsider]);

            do
            {
                nextMz = nextPeakMzs.Min();
                indexOfNextScanToConsider = Array.IndexOf(nextPeakMzs, nextMz);

                if ((nextMz - lastPeak.Mz) / lastPeak.Mz * 1e6 <= ppmTolerance)
                {
                    // Combine next peak with lastPeak
                    lastPeak.AddMzPeak(nextMz, nextPeaksIntensites[indexOfNextScanToConsider]);
                }
                else
                {
                    // Replace lastPeak
                    finalizedPeaks.Add(lastPeak);
                    lastPeak = new GeneratedPeak(nextMz, nextPeaksIntensites[indexOfNextScanToConsider]);
                }

                peaksLeft[indexOfNextScanToConsider]--;
                if (peaksLeft[indexOfNextScanToConsider] == 0)
                    nextPeakMzs[indexOfNextScanToConsider] = double.MaxValue;
                else
                {
                    nextPeakMzs[indexOfNextScanToConsider] = spectraToCombine[indexOfNextScanToConsider].XArray[totalPeaks[indexOfNextScanToConsider] - peaksLeft[indexOfNextScanToConsider]];
                    nextPeaksIntensites[indexOfNextScanToConsider] = spectraToCombine[indexOfNextScanToConsider].YArray[totalPeaks[indexOfNextScanToConsider] - peaksLeft[indexOfNextScanToConsider]];
                }
            } while (peaksLeft.Any(b => b > 0));

            finalizedPeaks.Add(lastPeak);

            return new GeneratedMzSpectrum(finalizedPeaks.Select(b => b.Mz).ToArray(), finalizedPeaks.Select(b => b.Intensity).ToArray(), false);
        }

        #endregion Private Methods
    }
}