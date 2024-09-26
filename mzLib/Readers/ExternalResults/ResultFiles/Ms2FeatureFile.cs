using CsvHelper;

namespace Readers
{
    /// <summary>
    /// Concrete Product for reading and representing a ms2.feature deconvolution result
    /// For supported versions and software this file type can come from see
    ///     Readers.ExternalResources.SupportedVersions.txt
    /// </summary>
    public class Ms2FeatureFile : ResultFile<Ms2Feature>, IResultFile
    {
        public override SupportedFileType FileType => SupportedFileType.Ms2Feature;

        public override Software Software { get; set; }

        public Ms2FeatureFile(string filePath, Software deconSoftware = Software.Unspecified) : base(filePath, deconSoftware) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public Ms2FeatureFile() : base() {}

        /// <summary>
        /// Load Results to the Results List from the given filepath
        /// </summary>
        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), Ms2Feature.CsvConfiguration);
            Results = csv.GetRecords<Ms2Feature>().ToList();
        }

        /// <summary>
        /// Writes results to a specific output path
        /// </summary>
        /// <param name="outputPath">destination path</param>
        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), Ms2Feature.CsvConfiguration);

            csv.WriteHeader<Ms2Feature>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }
    }
}
