using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;

namespace Readers.ReaderFactories
{
    /// <summary>
    /// Creates specific reader for MzML files
    /// </summary>
    internal sealed class MzMLReaderFactory : BaseReaderFactory, IReaderFactory
    {
        public MsDataFile Reader { get; }
        public MsDataFile CreateReader()
        {
            return new Mzml(FilePath);
        }

        internal MzMLReaderFactory(string filePath) : base(filePath)
        {
            FilePath = filePath;
            Reader = CreateReader();
        }
    }

    /// <summary>
    /// Creates specific reader for thermo raw files
    /// </summary>
    internal sealed class ThermoRawReaderFactory : BaseReaderFactory, IReaderFactory
    {
        public MsDataFile Reader { get; }

        internal ThermoRawReaderFactory(string filePath) : base(filePath)
        {
            FilePath = filePath;
            Reader = CreateReader();
        }
        public MsDataFile CreateReader()
        {
            return new ThermoRawFileReader(FilePath);
        }
    }

    /// <summary>
    /// Creates specific reader for Mgf files
    /// </summary>
    internal sealed class MgfReaderFactory : BaseReaderFactory, IReaderFactory
    {
        public MsDataFile Reader { get; }

        internal MgfReaderFactory(string filePath) : base(filePath)
        {
            FilePath = filePath;
            Reader = CreateReader();
        }

        public MsDataFile CreateReader()
        {
            return new Mgf(FilePath);
        }
    }
}
