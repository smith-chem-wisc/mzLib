using CsvHelper;
using System.IO;
using System.Linq;

namespace Readers;

public class FdrBenchPeptideFile : ResultFile<FdrBenchPeptide>
{
    public override SupportedFileType FileType => SupportedFileType.FdrBenchPeptide;
    public override Software Software { get; set; }

    public FdrBenchPeptideFile() : base()
    {
        Software = Software.FdrBench;
    }

    public FdrBenchPeptideFile(string filePath) : base(filePath, Software.FdrBench)
    {
    }

    public override void LoadResults()
    {
        using var reader = new StreamReader(FilePath);
        using var csv = new CsvReader(reader, FdrBenchPeptide.CsvConfiguration);
        Results = csv.GetRecords<FdrBenchPeptide>().ToList();
    }

    public override void WriteResults(string outputPath)
    {
        if (!CanRead(outputPath))
        {
            outputPath += FileType.GetFileExtension();
        }

        using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), FdrBenchPeptide.CsvConfiguration);
        csv.WriteHeader<FdrBenchPeptide>();
        foreach (var record in Results)
        {
            csv.NextRecord();
            csv.WriteRecord(record);
        }
    }
}
