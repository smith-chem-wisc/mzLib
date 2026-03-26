using System.Globalization;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;

namespace Readers;

public class FdrBenchProtein
{
    public static CsvConfiguration CsvConfiguration { get; } = new CsvConfiguration(CultureInfo.InvariantCulture)
    {
        Delimiter = "\t",
        HasHeaderRecord = true,
        TrimOptions = TrimOptions.Trim,
        BadDataFound = null,
    };

    [Name("protein")]
    public string Protein { get; set; } = string.Empty;

    [Name("q_value")]
    public double QValue { get; set; }

    [Name("score")]
    public double Score { get; set; }
}
