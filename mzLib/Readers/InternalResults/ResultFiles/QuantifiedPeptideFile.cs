using CsvHelper;
using MzLibUtil;
using System.Globalization;

namespace Readers;

/// <summary>
/// Reads a FlashLFQ peptide table (<c>QuantifiedPeptides.tsv</c>, or MetaMorpheus's
/// <c>AllQuantifiedPeptides.tsv</c>): the fixed columns by name and the per-sample <c>Intensity_</c>
/// and <c>Detection Type_</c> columns (and IsoTracker's <c>RetentionTime (min)_</c>) by header prefix.
/// </summary>
public class QuantifiedPeptideFile : ResultFile<QuantifiedPeptideFromTsv>, IResultFile
{
    private static readonly string[] SamplePrefixes = ["Intensity_", "Detection Type_", "RetentionTime (min)_"];

    public override SupportedFileType FileType => SupportedFileType.FlashLFQQuantifiedPeptide;
    public override Software Software { get; set; }

    public QuantifiedPeptideFile(string filePath) : base(filePath, Software.FlashLFQ) { }

    /// <summary>
    /// Constructor used to initialize from the factory method
    /// </summary>
    public QuantifiedPeptideFile() : base() { }

    /// <exception cref="MzLibException">The file is not a readable peptide table.</exception>
    public override void LoadResults()
    {
        try
        {
            using var csv = new CsvReader(new StreamReader(FilePath), QuantifiedPeptideFromTsv.CsvConfiguration);
            csv.Read();
            csv.ReadHeader();
            var header = csv.HeaderRecord ?? [];
            var sampleColumns = new List<(int Index, string Prefix, string Label)>();
            for (int i = 0; i < header.Length; i++)
            {
                string? prefix = SamplePrefixes.FirstOrDefault(p => header[i].StartsWith(p, StringComparison.Ordinal));
                if (prefix != null)
                    sampleColumns.Add((i, prefix, header[i][prefix.Length..]));
            }

            var results = new List<QuantifiedPeptideFromTsv>();
            while (csv.Read())
            {
                var row = csv.GetRecord<QuantifiedPeptideFromTsv>();
                var samples = new Dictionary<string, QuantifiedPeptideSample>();
                foreach (var (index, prefix, label) in sampleColumns)
                {
                    if (!samples.TryGetValue(label, out var s))
                        samples[label] = s = new QuantifiedPeptideSample(label);
                    string cell = csv.GetField(index) ?? "";
                    bool blank = string.IsNullOrWhiteSpace(cell);
                    switch (prefix)
                    {
                        case "Intensity_": s.Intensity = blank ? null : double.Parse(cell, NumberStyles.Float, CultureInfo.InvariantCulture); break;
                        case "RetentionTime (min)_": s.RetentionTime = blank ? null : double.Parse(cell, NumberStyles.Float, CultureInfo.InvariantCulture); break;
                        default: s.DetectionType = blank ? null : cell; break;
                    }
                }
                row.Samples = samples;
                results.Add(row);
            }
            Results = results;
        }
        catch (Exception e) when (e is not MzLibException)
        {
            throw new MzLibException($"Could not read peptide file '{FilePath}': {e.Message}", e);
        }
    }

    /// <exception cref="NotSupportedException">Always. FlashLFQ writes this format from its own peptides.</exception>
    public override void WriteResults(string outputPath) =>
        throw new NotSupportedException(
            "Writing peptide tables is not supported; FlashLFQ writes them from FlashLFQ.Peptide.");
}
