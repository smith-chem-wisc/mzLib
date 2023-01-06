using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using MzLibUtil;

namespace MzLibSpectralAveraging
{
    /// <summary>
    /// Manages the spectra files being averaged
    /// </summary>
    public static class SpectraFileProcessing
    {
        public static MsDataScan[] ProcessSpectra(List<MsDataScan> scans, MzLibSpectralAveragingOptions options)
        {
            switch (options.SpectraFileProcessingType)
            {
                case SpectraFileProcessingType.AverageAll:
                    return AverageAll(scans, options.SpectralAveragingOptions);

                case SpectraFileProcessingType.AverageEverynScans:
                    options.ScanOverlap = 0;
                    return AverageEverynScans(scans, options);

                case SpectraFileProcessingType.AverageEverynScansWithOverlap:
                    return AverageEverynScans(scans, options);

                case SpectraFileProcessingType.AverageDDAScans:
                    options.ScanOverlap = 0;
                    return AverageDDAScans(scans, options);

                case SpectraFileProcessingType.AverageDDAScansWithOverlap:
                    return AverageDDAScans(scans, options);


                    // TODO: Implement these two
                case SpectraFileProcessingType.MovingAverageMs1Scans:
                    return AverageMovingAverageMs1Scans(scans, options);

                case SpectraFileProcessingType.MovingAverageAllScans:
                    return AverageMovingAverageAllScans(scans, options);


                default: throw new MzLibException("Averaging spectra file processing type not yet implemented");
            }
        }

        private static MsDataScan[] AverageAll(List<MsDataScan> scans, SpectralAveragingOptions options)
        {
            // average spectrum
            MsDataScan representativeScan = scans.First();
            MzSpectrum averagedSpectrum = scans.CombineSpectra(options);

            // create output
            MsDataScan averagedScan = new(averagedSpectrum, 1, representativeScan.OneBasedScanNumber,
                representativeScan.IsCentroid, representativeScan.Polarity, scans.Select(p => p.RetentionTime).Average(),
                averagedSpectrum.Range, null, representativeScan.MzAnalyzer, scans.Select(p => p.TotalIonCurrent).Average(),
                representativeScan.InjectionTime, null, representativeScan.NativeId);
            MsDataScan[] msDataScans = new MsDataScan[] { averagedScan };
            return msDataScans;
        }

        private static MsDataScan[] AverageEverynScans(List<MsDataScan> scans, MzLibSpectralAveragingOptions options)
        {
            List<MsDataScan> averagedScans = new();
            int scanNumberIndex = 1;
            for (int i = 0; i < scans.Count; i += options.NumberOfScansToAverage - options.ScanOverlap)
            {
                // get the scans to be averaged
                List<MsDataScan> scansToProcess = new();
                if (i <= options.ScanOverlap) // very start of the file
                {
                    scansToProcess = scans.GetRange(i, options.NumberOfScansToAverage);
                }
                else if (i + options.NumberOfScansToAverage > scans.Count) // very end of the file
                {
                    break;
                }
                else // anywhere in the middle of the file
                {
                    scansToProcess = scans.GetRange(i , options.NumberOfScansToAverage);
                }

                // average scans
                MsDataScan representativeScan = scansToProcess.First();
                MzSpectrum averagedSpectrum = scansToProcess.CombineSpectra(options.SpectralAveragingOptions);
                MsDataScan averagedScan = new(averagedSpectrum, scanNumberIndex, 1,
                    representativeScan.IsCentroid, representativeScan.Polarity, scansToProcess.Select(p => p.RetentionTime).Minimum(),
                    averagedSpectrum.Range, null, representativeScan.MzAnalyzer, scansToProcess.Select(p => p.TotalIonCurrent).Average(),
                    scansToProcess.Select(p => p.InjectionTime).Average(), null, representativeScan.NativeId);
                string newNativeID =
                    averagedScan.NativeId.Replace(averagedScan.NativeId.Split("=").Last(), scanNumberIndex.ToString());
                averagedScan.SetNativeID(newNativeID);
                averagedScans.Add(averagedScan);
                scanNumberIndex++;
            }

            return averagedScans.ToArray();
        }

        private static MsDataScan[] AverageDDAScans(List<MsDataScan> scans, MzLibSpectralAveragingOptions options)
        {
            List<MsDataScan> averagedScans = new();
            List<MsDataScan> ms1Scans = scans.Where(p => p.MsnOrder == 1).ToList();
            List<MsDataScan> ms2Scans = scans.Where(p => p.MsnOrder == 2).ToList();
            List<MsDataScan> scansToProcess = new();

            int scanNumberIndex = 1;
            for (int i = 0; i < ms1Scans.Count; i += options.NumberOfScansToAverage - options.ScanOverlap)
            {
                // get the scans to be averaged
                scansToProcess.Clear();
                IEnumerable<MsDataScan> ms2ScansFromAveragedScans;
                if (i + options.NumberOfScansToAverage > ms1Scans.Count) // very end of the file
                {
                    break;
                }
                else 
                {
                    scansToProcess = ms1Scans.GetRange(i, options.NumberOfScansToAverage);
                    // if next iteration breaks the loop (end of file), then add the rest of the MS2's
                    if ((i + options.NumberOfScansToAverage - options.ScanOverlap) + options.NumberOfScansToAverage > ms1Scans.Count)
                        ms2ScansFromAveragedScans = ms2Scans.Where(p =>
                            scansToProcess.Any(m => m.OneBasedScanNumber == p.OneBasedPrecursorScanNumber));
                    // if not, add MS2 scans from MS1's that will not be averaged in the next iteration
                    else
                        ms2ScansFromAveragedScans = ms2Scans.Where(p =>
                            scansToProcess.GetRange(0, options.NumberOfScansToAverage - options.ScanOverlap).Any(m => m.OneBasedScanNumber == p.OneBasedPrecursorScanNumber));
                }

                // average scans and add to averaged list
                MsDataScan representativeScan = scansToProcess.First();
                MzSpectrum averagedSpectrum = scansToProcess.CombineSpectra(options.SpectralAveragingOptions);
                MsDataScan averagedScan = new(averagedSpectrum, scanNumberIndex, 1,
                    representativeScan.IsCentroid, representativeScan.Polarity, representativeScan.RetentionTime,
                    averagedSpectrum.Range, null, representativeScan.MzAnalyzer, scansToProcess.Select(p => p.TotalIonCurrent).Average(),
                    scansToProcess.Select(p => p.InjectionTime).Average(), representativeScan.NoiseData, 
                    representativeScan.NativeId, representativeScan.SelectedIonMZ, representativeScan.SelectedIonChargeStateGuess, 
                    representativeScan.SelectedIonIntensity, representativeScan.IsolationMz, representativeScan.IsolationWidth, 
                    representativeScan.DissociationType, representativeScan.OneBasedPrecursorScanNumber, representativeScan.SelectedIonMonoisotopicGuessIntensity);
                string newNativeID =
                    averagedScan.NativeId.Replace(averagedScan.NativeId.Split("=").Last(), scanNumberIndex.ToString());
                averagedScan.SetNativeID(newNativeID);
                averagedScans.Add(averagedScan);
                int precursorScanIndex = scanNumberIndex;
                scanNumberIndex++;

                foreach (var scan in ms2ScansFromAveragedScans)
                {
                    MsDataScan newScan = new(scan.MassSpectrum, scanNumberIndex, scan.MsnOrder, scan.IsCentroid,
                        scan.Polarity, scan.RetentionTime, scan.ScanWindowRange, scan.ScanFilter, scan.MzAnalyzer,
                        scan.TotalIonCurrent, scan.InjectionTime, scan.NoiseData, scan.NativeId, scan.SelectedIonMZ, 
                        scan.SelectedIonChargeStateGuess, scan.SelectedIonIntensity, scan.IsolationMz, scan.IsolationWidth, 
                        scan.DissociationType, precursorScanIndex, scan.SelectedIonMonoisotopicGuessMz);
                    newNativeID =
                        newScan.NativeId.Replace(newScan.NativeId.Split("=").Last(), scanNumberIndex.ToString());
                    newScan.SetNativeID(newNativeID);
                    averagedScans.Add(newScan);
                    scanNumberIndex++;
                }
            }
            
            return averagedScans.ToArray();
        }

        private static MsDataScan[] AverageMovingAverageAllScans(List<MsDataScan> scans, MzLibSpectralAveragingOptions options)
        {
            throw new NotImplementedException();
        }

        private static MsDataScan[] AverageMovingAverageMs1Scans(List<MsDataScan> scans, MzLibSpectralAveragingOptions options)
        {
            throw new NotImplementedException();
        }
    }
}
