namespace Readers.ReaderFactories
{
    internal interface IReaderFactory
    {
        public MsDataFile Reader { get; }
        public MsDataFile CreateReader();
    }
}
