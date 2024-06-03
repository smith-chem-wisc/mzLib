
namespace Readers
{
    public class CruxResultFile : ResultFile<CruxResult>, IResultFile
    {
        public override SupportedFileType FileType => SupportedFileType.CruxResult;
        public override Software Software { get; set; }

        public CruxResultFile(string filePath) : base(filePath, Software.Crux) { }

        public CruxResultFile() : base() { }

        public override void LoadResults()
        {
            using var csv = new CsvHelper.CsvReader(new StreamReader(FilePath), CruxResult.CsvConfiguration);
            Results = csv.GetRecords<CruxResult>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            if (!CanRead(FilePath))
                outputPath += FileType.GetFileExtension();

            using (var csv = new CsvHelper.CsvWriter(new StreamWriter(File.Create(outputPath)), CruxResult.CsvConfiguration))
            {
                csv.WriteHeader<CruxResult>();
                foreach (var result in Results)
                {
                    csv.NextRecord();
                    csv.WriteRecord(result);
                }
            }
        }
    }
}
