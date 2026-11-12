using System.Globalization;
using System.Text;
using CsvHelper;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;

namespace Readers.DeNovoMapping;

/// <summary>
/// The tab-separated mapping summary: one <see cref="DeNovoMappingResult"/> per database, followed by
/// one prediction-count column per mzTab file.
/// </summary>
public class DeNovoMappingResultFile : ResultFile<DeNovoMappingResult>, IResultFile
{
    public static CsvConfiguration CsvConfiguration => new(CultureInfo.InvariantCulture)
    {
        Encoding = Encoding.UTF8,
        HasHeaderRecord = true,
        Delimiter = "\t",
        HeaderValidated = null,
        MissingFieldFound = null,
    };

    /// <summary>
    /// The columns every mapping summary has, in the order they are written.
    /// </summary>
    private static readonly string[] FixedColumns = typeof(DeNovoMappingResult).GetProperties()
        .Where(p => !Attribute.IsDefined(p, typeof(IgnoreAttribute)))
        .Select(p => p.Name)
        .ToArray();

    public override SupportedFileType FileType => SupportedFileType.DenovoMappingResult;

    public override Software Software { get; set; }

    public DeNovoMappingResultFile(string filePath) : base(filePath, Software.Mzlib) { }

    /// <summary>
    /// Constructor used to initialize from the factory method
    /// </summary>
    public DeNovoMappingResultFile() : base()
    {
        // The base constructor leaves Software unspecified.
        Software = Software.Mzlib;
    }

    /// <summary>
    /// True when a header has the columns that identify a mapping summary, whatever the file is called.
    /// </summary>
    public static bool HasMappingSummaryColumns(IReadOnlyCollection<string> headers) =>
        headers.Contains(nameof(DeNovoMappingResult.DatabaseIdentifier))
        && headers.Contains(nameof(DeNovoMappingResult.TargetPredictionScores))
        && headers.Contains(nameof(DeNovoMappingResult.DecoyPredictionScores));

    public override void LoadResults()
    {
        using var csv = new CsvReader(new StreamReader(FilePath), CsvConfiguration);
        if (!csv.Read())
        {
            Results = new();
            return;
        }

        csv.ReadHeader();
        var headers = csv.HeaderRecord ?? [];

        // Any column that is not a fixed one holds per-file prediction counts. Older MetaCurator summaries
        // also carry entrapment columns; those are not counts per file.
        var fileColumns = headers
            .Where(h => !FixedColumns.Contains(h) && !h.Contains("Entrapment", StringComparison.Ordinal))
            .ToArray();

        var results = new List<DeNovoMappingResult>();
        while (csv.Read())
        {
            var result = csv.GetRecord<DeNovoMappingResult>();
            foreach (var column in fileColumns)
            {
                if (csv.TryGetField(column, out int count))
                    result.PredictionCountsByFile[column] = count;
            }
            results.Add(result);
        }

        Results = results;
    }

    /// <summary>
    /// Writes the summary, sorted by total then target predictions. A path that does not end in
    /// MappingSummary.tsv gets that appended.
    /// </summary>
    public override void WriteResults(string outputPath)
    {
        if (!CanRead(outputPath))
            outputPath += FileType.GetFileExtension();

        using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), CsvConfiguration);
        var fileColumns = Results
            .SelectMany(r => r.PredictionCountsByFile.Keys)
            .Distinct()
            .OrderBy(k => k, StringComparer.Ordinal)
            .ToList();

        csv.WriteHeader<DeNovoMappingResult>();
        foreach (var column in fileColumns)
            csv.WriteField(column);

        foreach (var result in Results.OrderByDescending(p => p.TotalPredictions).ThenByDescending(p => p.TargetPredictions))
        {
            result.FinalizeValues();

            csv.NextRecord();
            csv.WriteRecord(result);
            foreach (var column in fileColumns)
            {
                result.PredictionCountsByFile.TryGetValue(column, out int count);
                csv.WriteField(count);
            }
        }

        csv.NextRecord();
        FilePath = outputPath;
    }
}
