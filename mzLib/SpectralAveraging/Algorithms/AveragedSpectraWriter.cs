#nullable enable
using System;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using Readers;

namespace SpectralAveraging;
public static class AveragedSpectraWriter
{
    /// <summary>
    /// Public entry point for outputting averaged spectra
    /// </summary>
    /// <param name="averagedScans"></param>
    /// <param name="parameters"></param>
    /// <param name="originalSpectraPath"></param>
    /// <param name="destinationDirectory"></param>
    /// <param name="averagedFileName"></param>
    /// <remarks>
    ///     Destination directory and averaged file name can either both be set, one or the other, or none at all
    ///     If the destination directory and averaged file name are not specified,
    ///         the averaged spectra will be written to the same directory as the original
    ///         file path with the name $"{Original File Name}-averaged.mzML"
    ///     Destination directory can be set to allow averaged spectra to be sent to the desired location
    ///     Averaged file name will allow the mzML file to have a customized name
    /// </remarks>
    /// <exception cref="NotImplementedException"></exception>
    public static void WriteAveragedScans(MsDataScan[] averagedScans, SpectralAveragingParameters parameters,
        string originalSpectraPath, string? destinationDirectory = null, string? averagedFileName = null)
    {
        switch (parameters.OutputType)
        {
            case OutputType.MzML:
                WriteAveragedSpectraAsMzMl(averagedScans, originalSpectraPath, destinationDirectory, averagedFileName);
                break;

            default: throw new MzLibException("Output averaged scans type not implemented");
        }
    }

    /// <summary>
    /// Outputs averaged spectra in mzml format to same location as original spectra unless
    /// specified by the optional destination directory and destination file name fields
    /// </summary>
    /// <param name="averagedScans"></param>
    /// <param name="originalSpectraPath"></param>
    /// <param name="destinationDirectory"></param>
    /// <param name="averagedFileName"></param>
    private static void WriteAveragedSpectraAsMzMl(MsDataScan[] averagedScans,
        string originalSpectraPath, string? destinationDirectory = null, string? averagedFileName = null)
    {
        // gather necessary information to generate file from original spectra path
        var spectraDirectory = Path.GetDirectoryName(originalSpectraPath) ??
                               throw new MzLibException("Cannot Access Spectra Directory");
                               
        var sourceFile = MsDataFileReader.GetDataFile(originalSpectraPath).GetSourceFile();
        GenericMsDataFile msDataFile = new(averagedScans, sourceFile);

        // construct output path
        string averagedOutputPath;
        if (destinationDirectory is null)
        {
            // destination and file name not specified
            if (averagedFileName is null)
            {
                averagedOutputPath = Path.Combine(spectraDirectory,
                    PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(
                        originalSpectraPath) + "-averaged");
            }
            // destination not specified but name is
            else
            {
                averagedOutputPath = Path.Combine(spectraDirectory, averagedFileName);
            }
        }
        else
        {
            if (!Directory.Exists(destinationDirectory))
                Directory.CreateDirectory(destinationDirectory);

            // destination specified but not name
            if (averagedFileName is null)
            {
                averagedOutputPath = Path.Combine(destinationDirectory,
                    PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(
                        originalSpectraPath) + "-averaged");
            }
            // both destination and name specified
            else
            {
                averagedOutputPath = Path.Combine(destinationDirectory, averagedFileName);
            }
        }

        // add file extension
        if (!averagedOutputPath.EndsWith(".mzML"))
            averagedOutputPath += ".mzML";

        // check to see if file already exists, if so, add a number to the end
        int index = 1;
        while (File.Exists(averagedOutputPath))
        {
            int indexToInsert = averagedOutputPath.IndexOf(".mzML", StringComparison.InvariantCulture);
            
            // if first time needing to add an integer to filename
            averagedOutputPath = index == 1
                ? averagedOutputPath.Insert(indexToInsert, $"({index})")
                : averagedOutputPath.Replace($"({index - 1})", $"({index})");
            
            index++;
        }

        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(msDataFile, averagedOutputPath, true);
    }
}