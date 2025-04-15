using CsvHelper;
using Easy.Common.Extensions;

namespace Readers
{
    public class MsPathFinderTResultFile : ResultFile<MsPathFinderTResult>, IResultFile
    {
        public override SupportedFileType FileType => FilePath.ParseFileType();
        public override Software Software { get; set; }

        public MsPathFinderTResultFile(string filePath) : base(filePath, Software.MsPathFinderT) { }

        public MsPathFinderTResultFile() : base() { }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), MsPathFinderTResult.CsvConfiguration);
            Results = csv.GetRecords<MsPathFinderTResult>().ToList();
            if (Results.Any() && Results.First().FileNameWithoutExtension.IsNullOrEmpty())
                Results.ForEach(p => p.FileNameWithoutExtension = string.Join("_", Path.GetFileNameWithoutExtension(FilePath).Split('_')[..^1]));
        }

        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using (var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), MsPathFinderTResult.CsvConfiguration))
            {
                csv.WriteHeader<MsPathFinderTResult>();
                foreach (var result in Results)
                {
                    csv.NextRecord();
                    csv.WriteRecord(result);
                }
            }
        }
    }
}
