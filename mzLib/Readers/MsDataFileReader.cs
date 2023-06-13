using MassSpectrometry;
using MzLibUtil;

namespace Readers
{
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
                _ => throw new MzLibException("File extension not supported."),
            };
        }
    }
}
