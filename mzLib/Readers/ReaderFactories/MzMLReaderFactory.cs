namespace Readers.ReaderFactories
{
    internal class MzMLReaderFactory : BaseReaderFactory, IReaderFactory
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
}
