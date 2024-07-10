using CsvHelper;
using Easy.Common.Extensions;

namespace Readers
{ 
    public class MsFraggerExperimentFile : ResultFile<MsFraggerExperiment>, IResultFile
    {
        public override SupportedFileType FileType => SupportedFileType.MsFraggerExperiment;
        public override Software Software { get; set; }

        public MsFraggerExperimentFile(string filePath) : base(filePath, Software.Unspecified) { }

        /// <summary>
        ///        /// Constructor used to initialize from the factory method
        ///               /// </summary>
        ///                      public MsFraggerExperimentFile() : base() { }
        ///                      
        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), MsFraggerExperiment.CsvConfiguration);
            Results = csv.GetRecords<MsFraggerExperiment>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            throw new NotImplementedException();
        }
    }
}
