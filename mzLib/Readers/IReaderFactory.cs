using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;

namespace Readers
{
    public interface IReaderFactory
    {
        MsDataFile Reader { get; }
        MsDataFile CreateReader();
    }

    public class ReaderBase
    {
        public string FilePath { get; set; }
        public ReaderBase(string filePath)
        {
            FilePath = filePath;
        }
    }
    public class ThermoRawReaderFactory : ReaderBase, IReaderFactory
    {
        public MsDataFile Reader { get; private set; }

        public ThermoRawReaderFactory(string filePath) : base(filePath)
        {
            FilePath = filePath;
            Reader = CreateReader(); 
        } 
        public MsDataFile CreateReader()
        {
            return new ThermoRawFileReader(FilePath);
        }
    }

    public class MgfReaderFactory : ReaderBase, IReaderFactory
    {
        public MsDataFile Reader { get; }

        public MgfReaderFactory(string filePath) : base(filePath)
        {
            FilePath = filePath;
            Reader = CreateReader();
        }
        public MsDataFile CreateReader()
        {
            return new Mgf(FilePath); 
        }
    }

    public class MzMLReaderFactory : ReaderBase, IReaderFactory
    {
        public MsDataFile Reader { get; }
        public MsDataFile CreateReader()
        {
            return new Mzml(FilePath); 
        }

        public MzMLReaderFactory(string filePath) : base(filePath)
        {
            FilePath = filePath;
            Reader = CreateReader(); 
        }
    }

    public class BrukerFileReader : ReaderBase, IReaderFactory
    {
        public MsDataFile Reader { get; }
        public MsDataFile CreateReader()
        {
            throw new NotImplementedException();
        }

        public BrukerFileReader(string filePath) : base(filePath)
        {
            throw new NotSupportedException();
        }
    }

}
