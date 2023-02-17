using System;
using System.IO;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;

namespace SpectralAveraging
{
    public static class AveragedSpectraWriter
    {
        /// <summary>
        ///     Public entry point for outputting averaged spectra
        /// </summary>
        /// <param name="averagedScans"></param>
        /// <param name="parameters"></param>
        /// <param name="originalSpectraPath"></param>
        /// <exception cref="NotImplementedException"></exception>
        public static void WriteAveragedScans(MsDataScan[] averagedScans, SpectralAveragingParameters parameters,
            string originalSpectraPath)
        {
            switch (parameters.OutputType)
            {
                case OutputType.MzML:
                    WriteAveragedSpectraAsMzML(averagedScans, originalSpectraPath);
                    break;

                case OutputType.Text:
                    WriteAveragedSpectraAsTxtFile(averagedScans, parameters, originalSpectraPath);
                    break;

                default: throw new MzLibException("Output averaged scans type not implemented");
            }
        }

        /// <summary>
        ///     Outputs averaged spectra in mzml format to same location
        /// </summary>
        /// <param name="averagedScans"></param>
        /// <param name="originalSpectraPath"></param>
        private static void WriteAveragedSpectraAsMzML(MsDataScan[] averagedScans,
            string originalSpectraPath)
        {
            var spectraDirectory = Path.GetDirectoryName(originalSpectraPath) ??
                                   throw new MzLibException("Cannot Access Spectra Directory");
            var sourceFile = SpectraFileHandler.GetSourceFile(originalSpectraPath);
            MsDataFile msDataFile = new(averagedScans, sourceFile);
            var averagedPath = Path.Combine(spectraDirectory,
                "Averaged_" +
                PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(originalSpectraPath) +
                ".mzML");

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
            string originalSpectraPath)
        {
            var spectraDirectory = Path.GetDirectoryName(originalSpectraPath) ??
                                   throw new MzLibException("Cannot Access Spectra Directory");
            if (parameters.SpectraFileAveragingType != SpectraFileAveragingType.AverageAll)
            {
                spectraDirectory = Path.Combine(spectraDirectory, "AveragedSpectra");
                if (!Directory.Exists(spectraDirectory))
                    Directory.CreateDirectory(spectraDirectory);
            }


            var averagedPath = Path.Combine(spectraDirectory,
                "Averaged_" +
                PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(originalSpectraPath) +
                ".txt");

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
}