using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IO.BrukerFileReader;
using MassSpectrometry;
using MzLibUtil;

namespace Readers
{
    public static class MsDataFileReader
    {
        public static MsDataFile GetDataFile(string filePath)
        {
            string fileExtension = Path.GetExtension(filePath).ToLowerInvariant();
            return fileExtension switch
            {
                ".raw" => new ThermoRawFileReader(filePath),
                ".mzml" => new Mzml(filePath),
                ".mgf" => new Mgf(filePath),
                ".d" => new Bruker(filePath), 
                _ => throw new MzLibException("File extension not supported."),
            };
        }
    }
}
