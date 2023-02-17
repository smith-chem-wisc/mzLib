using System;
using System.Collections.Generic;
using IO.MzML;
using IO.ThermoRawFileReader;
using MassSpectrometry;
using MzLibUtil;

namespace SpectralAveraging
{
    public static class SpectraFileHandler
    {
        /// <summary>
        ///     Creates a List of MsDataScans from a spectra file
        ///     Currently supports MzML and raw
        /// </summary>
        /// <param name="filepath"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException">thrown if file type is not supported</exception>
        public static List<MsDataScan> LoadAllScansFromFile(string filepath)
        {
            List<MsDataScan> scans;
            if (filepath.EndsWith(".mzML", StringComparison.InvariantCultureIgnoreCase))
            {
                var temp = Mzml.LoadAllStaticData(filepath);
                scans = temp.GetAllScansList();
            }
            else if (filepath.EndsWith(".raw", StringComparison.InvariantCultureIgnoreCase))
            {
                var temp = ThermoRawFileReader.LoadAllStaticData(filepath);
                scans = temp.GetAllScansList();
            }
            else
            {
                throw new MzLibException("Cannot load spectra");
            }

            return scans;
        }

        /// <summary>
        ///     Gets the source file for the spectra file at designated path
        /// </summary>
        /// <param name="filepath"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException">throws if file type is not supported</exception>
        public static SourceFile GetSourceFile(string filepath)
        {
            if (filepath.EndsWith(".mzML", StringComparison.InvariantCultureIgnoreCase))
                return Mzml.LoadAllStaticData(filepath).SourceFile;
            if (filepath.EndsWith(".raw", StringComparison.InvariantCultureIgnoreCase))
                return ThermoRawFileReader.LoadAllStaticData(filepath).SourceFile;
            throw new MzLibException("Cannot access SourceFile");
        }
    }
}