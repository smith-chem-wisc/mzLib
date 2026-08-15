using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Readers;
using MassSpectrometry;
using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// Entry point for reading a mass spectrometry data file. Dispatches on the file's extension and
    /// returns the appropriate reader, so callers do not need to know whether they were handed an mzML,
    /// a Thermo .raw, an MGF or a Bruker .d.
    /// </summary>
    public static class MsDataFileReader 
    {

        public static MsDataFile GetDataFile(string filePath)
        {
            return filePath.ParseFileType() switch
            {
                SupportedFileType.ThermoRaw => new ThermoRawFileReader(filePath),
                SupportedFileType.MzML => new Mzml(filePath),
                SupportedFileType.Mgf => new Mgf(filePath),
                SupportedFileType.BrukerD => new BrukerFileReader(filePath), 
                SupportedFileType.BrukerTimsTof => new TimsTofFileReader(filePath),
                SupportedFileType.Ms1Align => new Ms1Align(filePath),
                SupportedFileType.Ms2Align => new Ms2Align(filePath),
                _ => throw new MzLibException("File type not supported"),
            };
        }
    }
}
