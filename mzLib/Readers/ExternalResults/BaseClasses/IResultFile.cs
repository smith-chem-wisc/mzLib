namespace Readers
{
    /// <summary>
    /// Abstract product interface for all result files that can be read through factory methods
    /// </summary>
    /// <typeparam name="TFactory"></typeparam>
    public interface IResultFile
    {
        public string FilePath { get; internal set; }
        public SupportedFileType FileType { get; }
        public Software Software { get; set; }

        /// <summary>
        /// Load Results to the Results List from the given filepath
        /// </summary>
        public void LoadResults();

        /// <summary>
        /// Writes results to a specific output path
        /// </summary>
        /// <param name="outputPath">destination path</param>
        public void WriteResults(string outputPath);


    }
}
