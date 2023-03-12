namespace Readers.ReaderFactories
{
    internal class BrukerReaderFactory : BaseReaderFactory, IReaderFactory
    {
	    private MsDataFile _reader;
	    internal BrukerReaderFactory(string filePath) : base(filePath)
        {
            FilePath = filePath;
            _reader = CreateReader();
        }

	    MsDataFile IReaderFactory.Reader => _reader;

	    public MsDataFile CreateReader()
	    {
		    return new Bruker.Bruker(FilePath); 
	    }
    }
}
