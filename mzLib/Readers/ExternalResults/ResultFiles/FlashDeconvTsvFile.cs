using CsvHelper;

namespace Readers
{
    /// <summary>
    /// Concrete Product for reading and representing a FlashDeconv's ms1.tsv output file
    /// For supported versions and software this file type can come from see
    ///     Readers.ExternalResources.SupportedVersions.txt
    /// </summary>
    public class FlashDeconvTsvFile : ResultFile<FlashDeconvTsv>, IResultFile
    {
        public override SupportedFileType FileType => SupportedFileType.Tsv_FlashDeconv;

        public override Software Software { get; set; } = Software.FLASHDeconv;

        public FlashDeconvTsvFile(string filePath) : base(filePath, Software.FLASHDeconv) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public FlashDeconvTsvFile() : base() { }

        /// <summary>
        /// Load Results to the Results List from the given filepath
        /// </summary>
        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), FlashDeconvTsv.CsvConfiguration);
            Results = csv.GetRecords<FlashDeconvTsv>().ToList();
        }

        /// <summary>
        /// Writes results to a specific output path
        /// </summary>
        /// <param name="outputPath">destination path</param>
        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), FlashDeconvTsv.CsvConfiguration);

            csv.WriteHeader<FlashDeconvTsv>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }
    }
}
