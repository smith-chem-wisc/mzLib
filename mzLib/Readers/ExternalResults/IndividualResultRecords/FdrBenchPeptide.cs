using System.Globalization;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;

namespace Readers;

public class FdrBenchPeptide
{
    public static CsvConfiguration CsvConfiguration { get; } = new CsvConfiguration(CultureInfo.InvariantCulture)
    {
        Delimiter = "\t",
        HasHeaderRecord = true,
        TrimOptions = TrimOptions.Trim,
        BadDataFound = null,
    };

    [Name("run")]
    public string Run { get; set; } = string.Empty;

    [Name("peptide")]
    public string Peptide { get; set; } = string.Empty;

    [Name("mod_peptide")]
    public string ModifiedPeptide { get; set; } = string.Empty;

    [Name("charge")]
    public int Charge { get; set; }

    [Name("q_value")]
    public double QValue { get; set; }

    [Name("PEP")]
    public double? Pep { get; set; }

    [Name("protein")]
    public string Protein { get; set; } = string.Empty;

    [Name("score")]
    public double Score { get; set; }
}
