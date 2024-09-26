using CsvHelper;

namespace Readers
{
    /// <summary>
    /// Concrete Product for reading and representing a FlashDeconv's ms1.tsv output file
    /// For supported versions and software this file type can come from see
    ///     Readers.ExternalResources.SupportedVersions.txt
    /// </summary>
    public class FlashDeconvMs1TsvFile : ResultFile<FlashDeconvMs1Tsv>, IResultFile
    {
        public override SupportedFileType FileType => SupportedFileType.Ms1Tsv_FlashDeconv;
        public override Software Software { get; set; }

        public FlashDeconvMs1TsvFile(string filePath) : base(filePath, Software.FLASHDeconv) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public FlashDeconvMs1TsvFile() : base() { }

        /// <summary>
        /// Load Results to the Results List from the given filepath
        /// </summary>
        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), FlashDeconvMs1Tsv.CsvConfiguration);
            Results = csv.GetRecords<FlashDeconvMs1Tsv>().ToList();
        }

        /// <summary>
        /// Writes results to a specific output path
        /// </summary>
        /// <param name="outputPath">destination path</param>
        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), FlashDeconvMs1Tsv.CsvConfiguration);

            csv.WriteHeader<FlashDeconvMs1Tsv>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }
    }
}
