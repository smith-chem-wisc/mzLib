using System.Globalization;
using System.Text;
using CsvHelper;
using CsvHelper.Configuration;

namespace Readers;

public class ChimerysPsmFile : ResultFile<ChimerysPsm>, IResultFile
{
    public static CsvConfiguration CsvConfiguration => new CsvConfiguration(CultureInfo.InvariantCulture)
    {
        Encoding = Encoding.UTF8,
        HasHeaderRecord = true,
        Delimiter = "\t",
        IgnoreBlankLines = true,
        TrimOptions = TrimOptions.Trim
    };

    public override SupportedFileType FileType => SupportedFileType.ChimerysPsm;
    public override Software Software { get; set; } = Software.Chimerys;

    public ChimerysPsmFile() : base() { }
    public ChimerysPsmFile(string path) : base(path) { }

    public override void LoadResults()
    {
        using var csv = new CsvReader(new StreamReader(FilePath), CsvConfiguration);
        Results = csv.GetRecords<ChimerysPsm>().ToList();
    }

    public override void WriteResults(string outputPath)
    {
        using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), CsvConfiguration);

        csv.WriteHeader<ChimerysPsm>();
        foreach (var result in Results)
        {
            csv.NextRecord();
            csv.WriteRecord(result);
        }
    }
}