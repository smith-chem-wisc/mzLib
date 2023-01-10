using System;
using System.IO;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;

namespace MzLibSpectralAveraging
{
    public static class AveragedSpectraOutputter
    {
        /// <summary>
        /// Public entry point for outputting averaged spectra
        /// </summary>
        /// <param name="averagedScans"></param>
        /// <param name="options"></param>
        /// <param name="originalSpectraPath"></param>
        /// <exception cref="NotImplementedException"></exception>
        public static void OutputAveragedScans(MsDataScan[] averagedScans, SpectralAveragingOptions options,
           string originalSpectraPath)
        {
            
            switch (options.OutputType)
            {
                case OutputType.mzML:
                    OutputAveragedSpectraAsMzML(averagedScans, originalSpectraPath);
                    break;

                case OutputType.txt:
                    OutputAveragedSpectraAsTxtFile(averagedScans, options, originalSpectraPath);
                    break;

                default: throw new MzLibException("Output averaged scans type not implemented");
            }

        }

        /// <summary>
        /// Outputs averaged spectra in mzml format to same location
        /// </summary>
        /// <param name="averagedScans"></param>
        /// <param name="originalSpectraPath"></param>
        private static void OutputAveragedSpectraAsMzML(MsDataScan[] averagedScans,
            string originalSpectraPath)
        {
            var spectraDirectory = Path.GetDirectoryName(originalSpectraPath) ?? throw new MzLibException("Cannot Access Spectra Directory");
            SourceFile sourceFile = SpectraFileHandler.GetSourceFile(originalSpectraPath);
            MsDataFile msDataFile = new(averagedScans, sourceFile);
            string averagedPath = Path.Combine(spectraDirectory,
                "Averaged_" + PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(originalSpectraPath) + ".mzML");

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(msDataFile, averagedPath, true);
        }

        /// <summary>
        /// Outputs averaged spectra in text format, with each spectra being 
        /// </summary>
        /// <param name="averagedScans"></param>
        /// <param name="options"></param>
        /// <param name="originalSpectraPath"></param>
        /// <exception cref="MzLibException"></exception>
        private static void OutputAveragedSpectraAsTxtFile(MsDataScan[] averagedScans, SpectralAveragingOptions options,
            string originalSpectraPath)
        {
            var spectraDirectory = Path.GetDirectoryName(originalSpectraPath) ?? throw new MzLibException("Cannot Access Spectra Directory");
            if (options.SpectraFileProcessingType != SpectraFileProcessingType.AverageAll)
            {
                spectraDirectory = Path.Combine(spectraDirectory, "AveragedSpectra");
                if (!Directory.Exists(spectraDirectory))
                    Directory.CreateDirectory(spectraDirectory);
            }

            
            string averagedPath = Path.Combine(spectraDirectory,
                "Averaged_" + PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(originalSpectraPath)  + ".txt");

            foreach (var scan in averagedScans)
            {
                if (options.SpectraFileProcessingType != SpectraFileProcessingType.AverageAll)
                    averagedPath = Path.Combine(spectraDirectory,
                        "Averaged_" +
                        PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(originalSpectraPath) +
                        "_" + scan.OneBasedScanNumber + ".txt");
                using StreamWriter writer = new StreamWriter(File.Create(averagedPath));

                for (int i = 0; i < scan.MassSpectrum.XArray.Length; i++)
                {
                    writer.WriteLine(scan.MassSpectrum.XArray[i] + "," + scan.MassSpectrum.YArray[i]);
                }
            }
        }
    }
}
