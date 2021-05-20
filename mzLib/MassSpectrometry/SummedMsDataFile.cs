using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    public class SummedMsDataFile : MsDataFile
    {
        private readonly MsDataFile raw;
        private readonly int numScansToAverage;
        private readonly double ppmToleranceForPeakCombination;

        public SummedMsDataFile(MsDataFile raw, int numScansToAverage, double ppmToleranceForPeakCombination) :
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

        public override List<MsDataScan> GetAllScansList()
        {
            List<MsDataScan> allScans = new List<MsDataScan>();
            for (int scanNumber = 1; scanNumber <= Scans.Length; scanNumber++)
            {
                allScans.Add(GetOneBasedScan(scanNumber));
            }
            return allScans;
        }

        public override MsDataScan GetOneBasedScan(int oneBasedScanNumber)
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

                MzSpectrum peaks = CombinePeaks(raw.GetAllScansList().Where(b => b.OneBasedScanNumber >= oneBasedScanNumber && b.OneBasedScanNumber <= oneBasedScanNumber + numScansToAverage - 1).Select(b => b.MassSpectrum).ToList(), ppmToleranceForPeakCombination);

                MzRange scanWindowRange = representative.ScanWindowRange;

                double totalIonCurrent = peaks.SumOfAllY;
                double injectionTime = double.NaN;
                double[,] noiseData = null;

                Scans[oneBasedScanNumber - 1] = new MsDataScan(peaks, oneBasedScanNumber, msnOrder, isCentroid, polarity, retentionTime, scanWindowRange, null, mzAnalyzer, totalIonCurrent, injectionTime, noiseData, "scan=" + oneBasedScanNumber);
            }
            return Scans[oneBasedScanNumber - 1];
        }

        private static MzSpectrum CombinePeaks(List<MzSpectrum> spectraToCombine, double ppmTolerance)
        {
            List<MzPeak> finalizedPeaks = new List<MzPeak>();

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
    }
}