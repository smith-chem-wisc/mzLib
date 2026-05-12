using CsvHelper;

namespace Readers
{
    /// <summary>
    /// Concrete Product for reading and representing a ms1.feature deconvolution result
    /// For supported versions and software this file type can come from see
    ///     Readers.ExternalResources.SupportedVersions.txt
    /// </summary>
    public class Ms1FeatureFile : ResultFile<Ms1Feature>, IResultFile, IMs1FeatureFile
    {
        public override SupportedFileType FileType => SupportedFileType.Ms1Feature;

        public sealed override Software Software { get; set; }

        public IEnumerable<ISingleChargeMs1Feature> GetMs1Features() => Results.SelectMany(r => r.GetSingleChargeFeatures());

        public Ms1FeatureFile(string filePath, Software deconSoftware = Software.Unspecified) : base(filePath,
            deconSoftware)
        {
            using (var sr = new StreamReader(filePath))
            {
                string firstLine = sr.ReadLine() ?? "";
                if (firstLine.Contains("\tApex_intensity\t") || firstLine.Contains("\tIntensity_Apex\t"))
                    Software = Software.TopFD;
                else
                    Software = Software.FLASHDeconv;
            }
        }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public Ms1FeatureFile() : base() { }

        /// <summary>
        /// Load Results to the Results List from the given filepath
        /// </summary>
        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), Ms1Feature.CsvConfiguration);
            Results = csv.GetRecords<Ms1Feature>().ToList();

            Software = Results.All(p => p.IntensityApex == null) ? Software.FLASHDeconv : Software.TopFD;
        }

        /// <summary>
        /// Writes results to a specific output path
        /// </summary>
        /// <param name="outputPath">destination path</param>
        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), Ms1Feature.CsvConfiguration);

            csv.WriteHeader<Ms1Feature>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }
    }
}
