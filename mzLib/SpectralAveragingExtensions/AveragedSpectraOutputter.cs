using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using System.IO;

namespace SpectralAveragingExtensions
{
    public static class AveragedSpectraOutputter
    {
        /// <summary>
        /// Public entry point for outputting averaged spectra
        /// </summary>
        /// <param name="averagedScans"></param>
        /// <param name="options"></param>
        /// <param name="spectraPath"></param>
        /// <exception cref="NotImplementedException"></exception>
        public static void OutputAveragedScans(MsDataScan[] averagedScans, MzLibSpectralAveragingOptions options,
           string spectraPath)
        {
            
            switch (options.OutputType)
            {
                case OutputType.mzML:
                    OutputAveragedSpectraAsMzML(averagedScans, spectraPath);
                    break;

                case OutputType.txt:
                    OutputAveragedSpectraAsTxtFile(averagedScans, options, spectraPath);
                    break;

                default: throw new MzLibException("Output averaged scans type not implemented");
            }

        }

        /// <summary>
        /// Outputs averaged spectra in mzml format to same location
        /// </summary>
        /// <param name="averagedScans"></param>
        /// <param name="spectraPath"></param>
        private static void OutputAveragedSpectraAsMzML(MsDataScan[] averagedScans,
            string spectraPath)
        {
            var spectraDirectory = Path.GetDirectoryName(spectraPath) ?? throw new MzLibException("Cannot Access Spectra Directory");
            SourceFile sourceFile = SpectraFileHandler.GetSourceFile(spectraPath);
            MsDataFile msDataFile = new(averagedScans, sourceFile);
            string averagedPath = Path.Combine(spectraDirectory,
                "Averaged_" + PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(spectraPath) + ".mzML");

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(msDataFile, averagedPath, true);
        }

        private static void OutputAveragedSpectraAsTxtFile(MsDataScan[] averagedScans, MzLibSpectralAveragingOptions options,
            string spectraPath)
        {
            var spectraDirectory = Path.GetDirectoryName(spectraPath) ?? throw new MzLibException("Cannot Access Spectra Directory");
            if (options.SpectraFileProcessingType != SpectraFileProcessingType.AverageAll)
            {
                spectraDirectory = Path.Combine(spectraDirectory, "AveragedSpectra");
                if (!Directory.Exists(spectraDirectory))
                    Directory.CreateDirectory(spectraDirectory);
            }

            
            string averagedPath = Path.Combine(spectraDirectory,
                "Averaged_" + PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(spectraPath)  + ".txt");

            foreach (var scan in averagedScans)
            {
                if (options.SpectraFileProcessingType != SpectraFileProcessingType.AverageAll)
                    averagedPath = Path.Combine(spectraDirectory,
                        "Averaged_" +
                        PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(spectraPath) +
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
