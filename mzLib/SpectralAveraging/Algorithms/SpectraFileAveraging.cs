using System;
using System.Collections.Generic;
using System.Linq;
using Easy.Common.Interfaces;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using MzLibUtil;

namespace SpectralAveraging;

/// <summary>
///     Manages the spectra files being averaged
/// </summary>
public static class SpectraFileAveraging
{
    /// <summary>
    ///     Averages an entire spectra file, grouping the spectra to be averaged based upon the
    ///     parameters.SpectraFileProcessing field
    /// </summary>
    /// <param name="scans">all scans from a file to be averaged</param>
    /// <param name="parameters">how to average the spectra</param>
    /// <returns></returns>
    /// <exception cref="MzLibException"></exception>
    public static MsDataScan[] AverageSpectraFile(List<MsDataScan> scans, SpectralAveragingParameters parameters)
    {
        switch (parameters.SpectraFileAveragingType)
        {
            case SpectraFileAveragingType.AverageAll:
                return AverageAll(scans, parameters);

            case SpectraFileAveragingType.AverageEverynScans:
                parameters.ScanOverlap = 0;
                return AverageEverynScans(scans, parameters);

            case SpectraFileAveragingType.AverageEverynScansWithOverlap:
                return AverageEverynScans(scans, parameters);

            case SpectraFileAveragingType.AverageDdaScans:
                return AverageDdaScans(scans, parameters);

            default: throw new MzLibException("Averaging spectra file processing type not yet implemented");
        }
    }

    /// <summary>
    ///     Average all scans
    /// </summary>
    /// <param name="scans">all scans from a file to be averaged</param>
    /// <param name="parameters">how to average the spectra</param>
    /// <returns></returns>
    private static MsDataScan[] AverageAll(IReadOnlyCollection<MsDataScan> scans, SpectralAveragingParameters parameters)
    {
        // average spectrum
        var representativeScan = scans.First();
        var averagedSpectrum = scans.AverageSpectra(parameters);

        // create output
        MsDataScan averagedScan = new(averagedSpectrum, 1, representativeScan.OneBasedScanNumber,
            representativeScan.IsCentroid, representativeScan.Polarity, scans.Select(p => p.RetentionTime).Average(),
            averagedSpectrum.Range, null, representativeScan.MzAnalyzer, scans.Select(p => p.TotalIonCurrent).Average(),
            representativeScan.InjectionTime, null, representativeScan.NativeId);
        MsDataScan[] msDataScans = { averagedScan };
        return msDataScans;
    }

    private static MsDataScan[] AverageEverynScans(List<MsDataScan> scans, SpectralAveragingParameters parameters)
    {
        List<MsDataScan> averagedScans = new();
        for (var i = 0; i < scans.Count; i += parameters.NumberOfScansToAverage - parameters.ScanOverlap)
        {
            // get the scans to be averaged
            List<MsDataScan> scansToProcess = new();
            if (i <= parameters.ScanOverlap) // very start of the file
                scansToProcess = scans.GetRange(i, parameters.NumberOfScansToAverage);
            else if (i + parameters.NumberOfScansToAverage > scans.Count) // very end of the file
                break;
            else // anywhere in the middle of the file
                scansToProcess = scans.GetRange(i, parameters.NumberOfScansToAverage);

            // average scans
            int middleIndex = scansToProcess.Count / 2;
            MsDataScan representativeScan = scansToProcess.Count % 2 == 0 ?
                scansToProcess[middleIndex - 1] :
                scansToProcess[middleIndex];

            var averagedSpectrum = scansToProcess.AverageSpectra(parameters);
            MsDataScan averagedScan = GetAveragedDataScanFromAveragedSpectrum(averagedSpectrum, representativeScan);
            averagedScans.Add(averagedScan);
        }

        return averagedScans.ToArray();
    }

    private static MsDataScan[] AverageDdaScans(List<MsDataScan> scans, SpectralAveragingParameters parameters)
    {
        List<MsDataScan> averagedScans = new();
        List<MsDataScan> scansToProcess = new();

        int representativeScanMs1Index = parameters.NumberOfScansToAverage % 2 == 0 ? // central scan
            parameters.NumberOfScansToAverage / 2 - 1 : parameters.NumberOfScansToAverage / 2;

        // iterate through all MS1 scans and average them
        foreach (var scan in scans.Where(p => p.MsnOrder == 1))
        {
            scansToProcess.Add(scan);
            // average with new scan from iteration, then remove first scan from list
            if (scansToProcess.Count != parameters.NumberOfScansToAverage) continue;

            MsDataScan centralScan = scansToProcess[representativeScanMs1Index];
            var averagedSpectrum = scansToProcess.AverageSpectra(parameters);
            var averagedScan = GetAveragedDataScanFromAveragedSpectrum(averagedSpectrum, centralScan);

            averagedScans.Add(averagedScan);
            scansToProcess.RemoveAt(0);
        }
        
        // add all scans that did not get averaged
        // this includes the MS1 scans from start and end of file and all MS2+ scans
        foreach (var unaveragedScan in scans.Where(original =>
                     !averagedScans.Select(avg => avg.OneBasedScanNumber).Contains(original.OneBasedScanNumber)))
            averagedScans.Add(unaveragedScan);
        
        return averagedScans.OrderBy(p => p.OneBasedScanNumber).ToArray();
    }

    private static MsDataScan GetAveragedDataScanFromAveragedSpectrum(MzSpectrum averagedSpectrum, 
        MsDataScan centralScan)
    {
        MsDataScan averagedScan = new(averagedSpectrum,
            centralScan.OneBasedScanNumber,
            1,
            centralScan.IsCentroid,
            centralScan.Polarity,
            centralScan.RetentionTime,
            averagedSpectrum.Range, null,
            centralScan.MzAnalyzer,
            averagedSpectrum.SumOfAllY,
            centralScan.InjectionTime,
            centralScan.NoiseData,
            centralScan.NativeId,
            centralScan.SelectedIonMZ,
            centralScan.SelectedIonChargeStateGuess,
            centralScan.SelectedIonIntensity,
            centralScan.IsolationMz,
            centralScan.IsolationWidth,
            centralScan.DissociationType,
            centralScan.OneBasedPrecursorScanNumber,
            centralScan.SelectedIonMonoisotopicGuessIntensity);
        var newNativeId =
            averagedScan.NativeId.Replace(averagedScan.NativeId.Split("=").Last(), centralScan.OneBasedScanNumber.ToString());
        averagedScan.SetNativeID(newNativeId);
        return averagedScan;
    }
}