namespace Readers.ReaderFactories
{
    internal class MgfReaderFactory : BaseReaderFactory, IReaderFactory
    {
        private MsDataFile _reader { get; }
        MsDataFile IReaderFactory.Reader => _reader; 

        internal MgfReaderFactory(string filePath) : base(filePath)
        {
            FilePath = filePath;
            _reader = CreateReader();
        }

        public MsDataFile CreateReader()
        {
            return new Mgf(FilePath);
        }
    }
}
