namespace Readers.ReaderFactories
{
    internal class MgfReaderFactory : BaseReaderFactory, IReaderFactory
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
