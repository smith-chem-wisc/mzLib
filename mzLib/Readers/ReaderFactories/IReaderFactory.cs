namespace Readers.ReaderFactories
{
    internal interface IReaderFactory
    {
        internal MsDataFile Reader { get; }
        public MsDataFile CreateReader();
    }
}
