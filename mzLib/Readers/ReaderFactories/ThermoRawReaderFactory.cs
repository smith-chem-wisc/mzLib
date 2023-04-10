namespace Readers.ReaderFactories
{
    internal class ThermoRawReaderFactory : BaseReaderFactory, IReaderFactory
    {
        public MsDataFile Reader { get; private set; }

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
}
