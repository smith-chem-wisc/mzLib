using System;
using System.Collections.Generic;
using IO.MzML;
using IO.ThermoRawFileReader;
using MassSpectrometry;
using MzLibUtil;

namespace MzLibSpectralAveraging
{
    public static class SpectraFileHandler
    {
        /// <summary>
        /// Creates a List of MsDataScans from a spectra file
        /// Currently supports MzML and raw
        /// </summary>
        /// <param name="filepath"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException">thrown if file type is not supported</exception>
        public static List<MsDataScan> LoadAllScansFromFile(string filepath)
        {
            List<MsDataScan> scans = new();
            if (filepath.EndsWith(".mzML"))
            {
                var temp = Mzml.LoadAllStaticData(filepath);
                scans = temp.GetAllScansList();
            }
            else if (filepath.EndsWith(".raw"))
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
        /// Gets the source file for the spectra file at designated path
        /// </summary>
        /// <param name="filepath"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException">throws if file type is not supported</exception>
        public static SourceFile GetSourceFile(string filepath)
        {
            List<MsDataScan> scans = new();
            if (filepath.EndsWith(".mzML"))
            {
                return Mzml.LoadAllStaticData(filepath).SourceFile;
            }
            else if (filepath.EndsWith(".raw"))
            {
                return ThermoRawFileReader.LoadAllStaticData(filepath).SourceFile;
            }
            else
            {
                throw new MzLibException("Cannot access SourceFile");
            }
        }
    }
}