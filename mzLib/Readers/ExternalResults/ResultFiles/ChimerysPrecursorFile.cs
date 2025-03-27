using System.Globalization;
using System.Text;
using CsvHelper;
using CsvHelper.Configuration;

namespace Readers;

public class ChimerysPrecursorFile : ResultFile<ChimerysPrecursor>, IResultFile
{
    public static CsvConfiguration CsvConfiguration => new CsvConfiguration(CultureInfo.InvariantCulture)
    {
        Encoding = Encoding.UTF8,
        HasHeaderRecord = true,
        Delimiter = "\t"
    };

    public override SupportedFileType FileType => SupportedFileType.ChimerysPrecursor;
    public override Software Software { get; set; } = Software.Chimerys;

    public ChimerysPrecursorFile() : base() { }
    public ChimerysPrecursorFile(string path) : base(path) { }

    public override void LoadResults()
    {
        using var csv = new CsvReader(new StreamReader(FilePath), CsvConfiguration);
        Results = csv.GetRecords<ChimerysPrecursor>().ToList();
    }

    public override void WriteResults(string outputPath)
    {
        if (!CanRead(outputPath))
            outputPath += FileType.GetFileExtension();

        using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), CsvConfiguration);

        csv.WriteHeader<ChimerysPrecursor>();
        foreach (var result in Results)
        {
            csv.NextRecord();
            csv.WriteRecord(result);
        }
    }
}