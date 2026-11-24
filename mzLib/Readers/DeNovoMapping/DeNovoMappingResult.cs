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

    #region Peptide level

    // One value per mapped target peptide rather than per prediction, so a peptide predicted many times counts once.
    // Always written; [Optional] so a summary written before they existed reads them as empty, or -1 for the count.

    /// <summary>
    /// For each target peptide with at least one prediction that has both retention times, the median of those
    /// predictions' <see cref="RetentionTimeErrors"/>. Peptides come in the order of their first prediction.
    /// </summary>
    [Optional]
    [TypeConverter(typeof(CommaDelimitedDoubleListConverter))]
    public List<double> PeptideRetentionTimeErrors { get; set; } = new();

    /// <summary>
    /// For each target peptide, the median score of its predictions, in the same peptide order.
    /// </summary>
    [Optional]
    [TypeConverter(typeof(CommaDelimitedDoubleListConverter))]
    public List<double> PeptidePredictionScores { get; set; } = new();

    /// <summary>
    /// Distinct decoy peptides matched by at least one prediction that no target in the database explains, or -1 when not
    /// known (a summary written before this column existed).
    /// </summary>
    [Optional]
    public int UniqueDecoyPeptidesMapped { get; set; } = -1;

    #endregion

    #region Host ambiguity

    // Filled only when the mapper was given host proteins (DeNovoMappingParameters.HostProteins). A hit is
    // host-ambiguous when the host also produces its peptide (DeNovoPrediction.IsHostShared) and unambiguous otherwise.
    // The fields above keep counting every hit. These columns ([HostColumn]) are written only when a result has a host,
    // and a summary without them reads them as -1, zero or empty.

    /// <summary>
    /// Target proteins in the host, or -1 when no host was given. Zero means a host was given but was empty, so every hit
    /// is unambiguous.
    /// </summary>
    [Optional, HostColumn]
    public int HostProteinCount { get; set; } = -1;

    /// <summary>
    /// Target predictions whose peptide the host does not produce. With <see cref="AmbiguousTargetPredictions"/>, sums to
    /// <see cref="TargetPredictions"/> when a host was given.
    /// </summary>
    [Optional, HostColumn]
    public int UnambiguousTargetPredictions { get; set; }

    /// <summary>
    /// Target predictions whose peptide the host also produces.
    /// </summary>
    [Optional, HostColumn]
    public int AmbiguousTargetPredictions { get; set; }

    /// <summary>
    /// Decoy predictions whose peptide the host does not produce, by the same rule as the targets, so unambiguous targets
    /// and decoys come from the same predictions.
    /// </summary>
    [Optional, HostColumn]
    public int UnambiguousDecoyPredictions { get; set; }

    /// <summary>
    /// Distinct target peptides matched by at least one unambiguous prediction.
    /// </summary>
    [Optional, HostColumn]
    public int UniqueUnambiguousPeptidesMapped { get; set; }

    /// <summary>
    /// Target proteins containing at least one peptide matched by an unambiguous prediction.
    /// </summary>
    [Optional, HostColumn]
    public int UniqueUnambiguousProteinsMapped { get; set; }

    /// <summary>
    /// <see cref="RetentionTimeErrors"/> of the unambiguous target predictions only.
    /// </summary>
    [Optional, HostColumn]
    [TypeConverter(typeof(CommaDelimitedDoubleListConverter))]
    public List<double> UnambiguousRetentionTimeErrors { get; set; } = new();

    [Optional, HostColumn]
    [TypeConverter(typeof(CommaDelimitedDoubleListConverter))]
    public List<double> UnambiguousTargetPredictionScores { get; set; } = new();

    [Optional, HostColumn]
    [TypeConverter(typeof(CommaDelimitedDoubleListConverter))]
    public List<double> UnambiguousDecoyPredictionScores { get; set; } = new();

    /// <summary>
    /// Distinct decoy peptides counted in <see cref="UniqueDecoyPeptidesMapped"/> that the host does not produce, or -1
    /// without a host.
    /// </summary>
    [Optional, HostColumn]
    public int UniqueUnambiguousDecoyPeptidesMapped { get; set; } = -1;

    /// <summary>
    /// <see cref="PeptideRetentionTimeErrors"/> of the peptides the host does not produce.
    /// </summary>
    [Optional, HostColumn]
    [TypeConverter(typeof(CommaDelimitedDoubleListConverter))]
    public List<double> UnambiguousPeptideRetentionTimeErrors { get; set; } = new();

    /// <summary>
    /// <see cref="PeptidePredictionScores"/> of the peptides the host does not produce.
    /// </summary>
    [Optional, HostColumn]
    [TypeConverter(typeof(CommaDelimitedDoubleListConverter))]
    public List<double> UnambiguousPeptidePredictionScores { get; set; } = new();

    /// <summary>
    /// True when the mapper was given host proteins, so the host-ambiguity fields are filled.
    /// </summary>
    [Ignore]
    public bool HasHost => HostProteinCount >= 0;

    #endregion

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

    internal static double Median(List<double> values)
    {
        var sorted = values.OrderBy(x => x).ToArray();
        int middle = sorted.Length / 2;
        return sorted.Length % 2 == 1 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
    }
}

/// <summary>
/// Marks a <see cref="DeNovoMappingResult"/> column that is written only when a result has a host.
/// </summary>
[AttributeUsage(AttributeTargets.Property)]
internal sealed class HostColumnAttribute : Attribute;

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