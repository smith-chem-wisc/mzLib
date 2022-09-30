using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;

namespace Readers
{
    public static class ReaderCreator
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
                    throw new NotImplementedException();
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
