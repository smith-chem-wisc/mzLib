namespace Readers.ReaderFactories
{
    internal class BaseReaderFactory
    {
        internal string FilePath { get; set; }
        internal BaseReaderFactory(string filePath)
        {
            FilePath = filePath;
        }
    }
}
