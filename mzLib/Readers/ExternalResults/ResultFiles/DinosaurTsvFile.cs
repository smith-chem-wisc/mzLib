using CsvHelper;

namespace Readers
{
    /// <summary>
    /// Concrete Product for reading and representing a Dinosaur .feature.tsv output file
    /// For supported versions and software this file type can come from see
    ///     Readers.ExternalResources.SupportedVersions.txt
    /// </summary>
    public class DinosaurTsvFile : ResultFile<DinosaurFeature>, IResultFile, IMs1FeatureFile
    {
        public override SupportedFileType FileType => SupportedFileType.Tsv_Dinosaur;
        public override Software Software { get; set; }
        public IEnumerable<ISingleChargeMs1Feature> GetMs1Features() { return Results; }

        public DinosaurTsvFile(string filePath) : base(filePath, Software.Dinosaur) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public DinosaurTsvFile() : base() { }

        /// <summary>
        /// Load Results to the Results List from the given filepath
        /// </summary>
        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), DinosaurFeature.CsvConfiguration);
            Results = csv.GetRecords<DinosaurFeature>().ToList();
        }

        /// <summary>
        /// Writes results to a specific output path
        /// </summary>
        /// <param name="outputPath">destination path</param>
        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), DinosaurFeature.CsvConfiguration);

            csv.WriteHeader<DinosaurFeature>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }
    }
}