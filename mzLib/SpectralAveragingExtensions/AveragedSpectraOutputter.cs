using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IO.MzML;
using MassSpectrometry;
using Nett;
using SpectralAveraging;
using SpectralAveragingExtensions.Util;

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

                default: throw new NotImplementedException("Output type not implemented");
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
            SourceFile sourceFile = SpectraFileHandler.GetSourceFile(spectraPath);
            MsDataFile msDataFile = new(averagedScans, sourceFile);
            string spectraDirectory = Path.GetDirectoryName(spectraPath);
            string averagedPath = Path.Combine(spectraDirectory,
                "Averaged_" + Path.GetFileNameWithoutExtension(spectraPath) + ".mzML");

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(msDataFile, averagedPath, true);
        }

        private static void OutputAveragedSpectraAsTxtFile(MsDataScan[] averagedScans, MzLibSpectralAveragingOptions options,
            string spectraPath)
        {
            string spectraDirectory = Path.GetDirectoryName(spectraPath);
            if (options.SpectraFileProcessingType != SpectraFileProcessingType.AverageAll)
            {
                spectraDirectory = Path.Combine(spectraDirectory, "AveragedSpectra");
                if (!Directory.Exists(spectraDirectory))
                    Directory.CreateDirectory(spectraDirectory);
            }

            string averagedPath = Path.Combine(spectraDirectory,
                "Averaged_" + Path.GetFileNameWithoutExtension(spectraPath)  + ".txt");

            foreach (var scan in averagedScans)
            {
                if (options.SpectraFileProcessingType != SpectraFileProcessingType.AverageAll)
                    averagedPath = Path.Combine(spectraDirectory,
                        "Averaged_" + Path.GetFileNameWithoutExtension(spectraPath) + "_" + scan.OneBasedScanNumber + ".txt");
                using (StreamWriter writer = new StreamWriter(File.Create(averagedPath)))
                {

                    for (int i = 0; i < scan.MassSpectrum.XArray.Length; i++)
                    {
                        writer.WriteLine(scan.MassSpectrum.XArray[i] + "," + scan.MassSpectrum.YArray[i]);
                    }
                }
            }
        }
    }
}
