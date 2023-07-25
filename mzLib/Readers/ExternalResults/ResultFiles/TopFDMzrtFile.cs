using CsvHelper;

namespace Readers
{
    /// <summary>
    /// Concrete Product for reading and representing a TopFD's mzrt.csv output file
    /// For supported versions and software this file type can come from see
    ///     Readers.ExternalResources.SupportedVersions.txt
    /// </summary>
    public class TopFDMzrtFile : ResultFile<TopFdMzrt>, IResultFile
    {
        public override SupportedFileType FileType => SupportedFileType.Mzrt_TopFd;

        public override Software Software { get; set; }

        public TopFDMzrtFile(string filePath) : base(filePath, Software.TopFD) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public TopFDMzrtFile() : base() { }

        /// <summary>
        /// Load Results to the Results List from the given filepath
        /// </summary>
        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), TopFdMzrt.CsvConfiguration);
            Results = csv.GetRecords<TopFdMzrt>().ToList();
        }

        /// <summary>
        /// Writes results to a specific output path
        /// </summary>
        /// <param name="outputPath">destination path</param>
        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), TopFdMzrt.CsvConfiguration);

            csv.WriteHeader<TopFdMzrt>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }
    }
}
