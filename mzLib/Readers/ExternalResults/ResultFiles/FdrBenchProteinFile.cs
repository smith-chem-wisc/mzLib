using CsvHelper;
using System.IO;

namespace Readers;

public class FdrBenchProteinFile : ResultFile<FdrBenchProtein>
{
    public override SupportedFileType FileType => SupportedFileType.FdrBenchProtein;
    public override Software Software { get; set; }

    public FdrBenchProteinFile() : base()
    {
        Software = Software.FdrBench;
    }

    public FdrBenchProteinFile(string filePath) : base(filePath, Software.FdrBench)
    {
    }

    public override void LoadResults()
    {
        throw new NotSupportedException("FDRBench protein exports are write-only.");
    }

    public override void WriteResults(string outputPath)
    {
        if (!CanRead(outputPath))
        {
            outputPath += FileType.GetFileExtension();
        }

        using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), FdrBenchProtein.CsvConfiguration);
        csv.WriteHeader<FdrBenchProtein>();
        foreach (var record in Results)
        {
            csv.NextRecord();
            csv.WriteRecord(record);
        }
    }
}
