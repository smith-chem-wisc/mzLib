namespace Readers.ReaderFactories
{
    internal class BrukerReaderFactory : BaseReaderFactory, IReaderFactory
    {
        public MsDataFile Reader { get; }

        internal BrukerReaderFactory(string filePath) : base(filePath)
        {
            FilePath = filePath;
            Reader = CreateReader();
        }
        public MsDataFile CreateReader()
        {
            throw new NotImplementedException();
        }
    }
}
