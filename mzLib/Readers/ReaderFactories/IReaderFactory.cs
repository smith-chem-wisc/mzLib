namespace Readers.ReaderFactories
{
    internal interface IReaderFactory
    {
        internal MsDataFile Reader { get; }
        internal MsDataFile CreateReader();
    }
}
