namespace Readers
{
    public class FileReader
    {
        public static TResultFile ReadFile<TResultFile>(string filePath) where TResultFile : IResultFile, new()
        {
            var resultFile = new TResultFile() { FilePath = filePath };
            resultFile.LoadResults();
            return resultFile;
        }
    }
}
