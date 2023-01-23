using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;
using Readers.ReaderFactories;

namespace Readers
{
    public static class MsDataFileReader
    {
        public static MsDataFile CreateReader(string filePath)
        {
            string fileExtension = Path.GetExtension(filePath).ToLowerInvariant();
            IReaderFactory factory = null; 
            switch (fileExtension)
            {
                case ".raw":
                    factory = new ThermoRawReaderFactory(filePath); 
                    break;
                case ".d":
                    factory = new BrukerReaderFactory(filePath);
                    break;
                case ".mzml":
                    factory = new MzMLReaderFactory(filePath); 
                    break;
                case ".mgf":
                    factory = new MgfReaderFactory(filePath); 
                    break;
                default:
                    throw new MzLibException("File extension not supported."); 
            }
            return factory.Reader; 
        }
    }
}
