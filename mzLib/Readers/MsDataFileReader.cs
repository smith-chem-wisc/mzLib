using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using MzLibUtil;
using Readers.ReaderFactories;

namespace Readers
{
    public static class MsDataFileReader
    {
        public static MsDataFile GetDataFile(string filePath)
        {
            string fileExtension = Path.GetExtension(filePath).ToLowerInvariant();
            IReaderFactory factory = null;
            factory = fileExtension switch
            {
                ".raw" => new ThermoRawReaderFactory(filePath),
                ".mzml" => new MzMLReaderFactory(filePath),
                ".mgf" => new MgfReaderFactory(filePath),
                _ => throw new MzLibException("File extension not supported."),
            };
            return factory.Reader; 
        }
    }
}
