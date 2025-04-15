using CsvHelper;

namespace Readers;

public class QuantifiedPeakFile : ResultFile<QuantifiedPeak>, IResultFile
{
    public override SupportedFileType FileType => SupportedFileType.FlashLFQQuantifiedPeak;
    public override Software Software { get; set; }

    public QuantifiedPeakFile(string filePath) : base(filePath, Software.FlashLFQ) { }

    /// <summary>
    /// Constructor used to initialize from the factory method
    /// </summary>
    public QuantifiedPeakFile() : base() { }

    public override void LoadResults()
    {
        using var csv = new CsvReader(new StreamReader(FilePath), QuantifiedPeak.CsvConfiguration);
        Results = csv.GetRecords<QuantifiedPeak>().ToList();
    }

    public override void WriteResults(string outputPath)
    {
        if (!CanRead(outputPath))
            outputPath += FileType.GetFileExtension();

        using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), QuantifiedPeak.CsvConfiguration);

        csv.WriteHeader<QuantifiedPeak>();
        foreach (var result in Results)
        {
            csv.NextRecord();
            csv.WriteRecord(result);
        }
    }
}