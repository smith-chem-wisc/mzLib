using CsvHelper;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using CsvHelper.TypeConversion;
using System.Globalization;

namespace Readers.DeNovoMapping;

/// <summary>
/// How the de novo predictions map onto one database: counts, score and retention-time distributions.
/// One row of the mapping summary file.
/// </summary>
public class DeNovoMappingResult
{
    public string DatabaseIdentifier { get; set; } = string.Empty;

    /// <summary>
    /// Target proteins in the database.
    /// </summary>
    public int DatabaseProteinCount { get; set; } = -1;

    /// <summary>
    /// Distinct target peptides produced by digesting the database, after I/L normalization.
    /// </summary>
    public int DatabasePeptideCount { get; set; } = -1;

    /// <summary>
    /// Target plus decoy predictions.
    /// </summary>
    public int TotalPredictions { get; set; }

    /// <summary>
    /// Predictions whose peptide is found in a target protein.
    /// </summary>
    public int TargetPredictions { get; set; }

    /// <summary>
    /// Predictions whose peptide is found in a decoy protein.
    /// </summary>
    public int DecoyPredictions { get; set; }

    /// <summary>
    /// Distinct target peptides matched by at least one prediction.
    /// </summary>
    public int UniquePeptidesMapped { get; set; }

    /// <summary>
    /// Target proteins containing at least one matched peptide.
    /// </summary>
    public int UniqueProteinsMapped { get; set; }

    public double MeanRtError { get; set; } = double.NaN;
    public double MedianRtError { get; set; } = double.NaN;
    public double StdDevRtError { get; set; } = double.NaN;
    public double MeanPredictionScore { get; set; } = double.NaN;
    public double MedianPredictionScore { get; set; } = double.NaN;
    public double StdDevPredictionScore { get; set; } = double.NaN;

    /// <summary>
    /// Observed minus predicted retention time, in minutes, for each target prediction that has both. The prediction is
    /// Chronologer's, mapped onto minutes by its file's calibration; predictions from files that could not be calibrated
    /// are left out.
    /// </summary>
    [TypeConverter(typeof(CommaDelimitedDoubleListConverter))]
    public List<double> RetentionTimeErrors { get; set; } = new();

    [TypeConverter(typeof(CommaDelimitedDoubleListConverter))]
    public List<double> PredictionScores { get; set; } = new();

    [TypeConverter(typeof(CommaDelimitedDoubleListConverter))]
    public List<double> TargetPredictionScores { get; set; } = new();

    [TypeConverter(typeof(CommaDelimitedDoubleListConverter))]
    public List<double> DecoyPredictionScores { get; set; } = new();

    /// <summary>
    /// Target plus decoy predictions per mzTab file. Written as one extra column per file.
    /// </summary>
    [Ignore]
    public Dictionary<string, int> PredictionCountsByFile { get; set; } = new();

    /// <summary>
    /// Recomputes the totals and summary statistics from the counts and distributions.
    /// </summary>
    public void FinalizeValues()
    {
        if (RetentionTimeErrors.Count > 0)
            (MeanRtError, MedianRtError, StdDevRtError) = Summarize(RetentionTimeErrors);

        TotalPredictions = TargetPredictions + DecoyPredictions;
        PredictionScores = TargetPredictionScores.Concat(DecoyPredictionScores).ToList();

        if (PredictionScores.Count > 0)
            (MeanPredictionScore, MedianPredictionScore, StdDevPredictionScore) = Summarize(PredictionScores);
    }

    private static (double Mean, double Median, double StdDev) Summarize(List<double> values)
    {
        double mean = values.Average();
        double median = Median(values);
        double stdDev = Math.Sqrt(values.Select(x => Math.Pow(x - mean, 2)).Average());
        return (mean, median, stdDev);
    }

    private static double Median(List<double> values)
    {
        var sorted = values.OrderBy(x => x).ToArray();
        int middle = sorted.Length / 2;
        return sorted.Length % 2 == 1 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
    }
}

/// <summary>
/// Reads and writes a list of doubles as one comma-separated field, culture-invariant.
/// </summary>
internal class CommaDelimitedDoubleListConverter : DefaultTypeConverter
{
    public override object ConvertFromString(string? text, IReaderRow row, MemberMapData memberMapData)
    {
        if (string.IsNullOrWhiteSpace(text))
            return new List<double>();

        return text.Split(',', StringSplitOptions.RemoveEmptyEntries)
            .Select(p => double.Parse(p, CultureInfo.InvariantCulture))
            .ToList();
    }

    public override string ConvertToString(object? value, IWriterRow row, MemberMapData memberMapData)
    {
        var values = value as IEnumerable<double> ?? [];
        return string.Join(',', values.Select(p => p.ToString("0.####", CultureInfo.InvariantCulture)));
    }
}