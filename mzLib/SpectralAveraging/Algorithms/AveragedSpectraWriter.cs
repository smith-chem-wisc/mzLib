using System;
using System.IO;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;

namespace SpectralAveraging;
public static class AveragedSpectraWriter
{
    /// <summary>
    ///     Public entry point for outputting averaged spectra
    /// </summary>
    /// <param name="averagedScans"></param>
    /// <param name="parameters"></param>
    /// <param name="originalSpectraPath"></param>
    /// <param name="destinationPath"></param>
    /// <exception cref="NotImplementedException"></exception>
    public static void WriteAveragedScans(MsDataScan[] averagedScans, SpectralAveragingParameters parameters,
        string originalSpectraPath, string? destinationPath = null)
    {
        switch (parameters.OutputType)
        {
            case OutputType.MzML:
                WriteAveragedSpectraAsMzMl(averagedScans, originalSpectraPath, destinationPath);
                break;

            case OutputType.Text:
                WriteAveragedSpectraAsTxtFile(averagedScans, parameters, originalSpectraPath, destinationPath);
                break;

            default: throw new MzLibException("Output averaged scans type not implemented");
        }
    }

    /// <summary>
    ///     Outputs averaged spectra in mzml format to same location
    /// </summary>
    /// <param name="averagedScans"></param>
    /// <param name="originalSpectraPath"></param>
    /// <param name="destinationPath"></param>
    private static void WriteAveragedSpectraAsMzMl(MsDataScan[] averagedScans,
        string originalSpectraPath, string? destinationPath = null)
    {
        var spectraDirectory = Path.GetDirectoryName(originalSpectraPath) ??
                               throw new MzLibException("Cannot Access Spectra Directory");
        var sourceFile = SpectraFileHandler.GetSourceFile(originalSpectraPath);
        MsDataFile msDataFile = new(averagedScans, sourceFile);

        string averagedPath;
        if (destinationPath == null)
            averagedPath = Path.Combine(spectraDirectory,
            "Averaged_" +
            PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(originalSpectraPath) +
            ".mzML");
        else
        {
            if (!destinationPath.EndsWith(".mzML"))
                averagedPath = destinationPath + ".mzML";
            else
                averagedPath = destinationPath;
        }

        int index = 1;
        while (File.Exists(averagedPath))
        {
            int indexToInsert = averagedPath.IndexOf(".mzML", StringComparison.InvariantCulture);
            averagedPath = averagedPath.Insert(indexToInsert, $"({index})");
        }

        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(msDataFile, averagedPath, true);
    }

    /// <summary>
    ///     Outputs averaged spectra in text format, with each spectra being
    /// </summary>
    /// <param name="averagedScans"></param>
    /// <param name="parameters"></param>
    /// <param name="originalSpectraPath"></param>
    /// <exception cref="MzLibException"></exception>
    private static void WriteAveragedSpectraAsTxtFile(MsDataScan[] averagedScans,
        SpectralAveragingParameters parameters,
        string originalSpectraPath, string? destinationPath = null)
    {
        var spectraDirectory = Path.GetDirectoryName(originalSpectraPath) ??
                               throw new MzLibException("Cannot Access Spectra Directory");
        if (parameters.SpectraFileAveragingType != SpectraFileAveragingType.AverageAll)
        {
            spectraDirectory = Path.Combine(spectraDirectory, "AveragedSpectra");
            if (!Directory.Exists(spectraDirectory))
                Directory.CreateDirectory(spectraDirectory);
        }

        string averagedPath;
        if (destinationPath == null)
            averagedPath = Path.Combine(spectraDirectory,
                "Averaged_" +
                PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(originalSpectraPath) +
                ".mzML");
        else
        {
            if (!destinationPath.EndsWith(".mzML"))
                averagedPath = destinationPath + ".mzML";
            else
                averagedPath = destinationPath;
        }

        foreach (var scan in averagedScans)
        {
            if (parameters.SpectraFileAveragingType != SpectraFileAveragingType.AverageAll)
                averagedPath = Path.Combine(spectraDirectory,
                    "Averaged_" +
                    PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(
                        originalSpectraPath) +
                    "_" + scan.OneBasedScanNumber + ".txt");
            using var writer = new StreamWriter(File.Create(averagedPath));

            for (var i = 0; i < scan.MassSpectrum.XArray.Length; i++)
                writer.WriteLine(scan.MassSpectrum.XArray[i] + "," + scan.MassSpectrum.YArray[i]);
        }
    }
}