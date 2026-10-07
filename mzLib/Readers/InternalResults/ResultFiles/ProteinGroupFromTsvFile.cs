using CsvHelper;
using MzLibUtil;
using System.Globalization;

namespace Readers;

/// <summary>
/// Reads a MetaMorpheus protein-group table, <c>AllQuantifiedProteinGroups.tsv</c>: the fixed columns
/// by name, and the per-sample-group <c>SpectralCount_</c>, <c>Intensity_</c>, <c>CountOccupancy_</c>
/// and <c>IntensityOccupancy_</c> columns by header prefix into
/// <see cref="ProteinGroupFromTsv.SampleGroups"/>.
/// </summary>
/// <remarks>
/// A sample group's block is two to four columns wide, and an isobaric file interleaves one counting
/// pair with one intensity pair per channel, so the columns are matched by prefix, never by position.
/// </remarks>
public class ProteinGroupFromTsvFile : ResultFile<ProteinGroupFromTsv>, IResultFile
{
    // No prefix is a prefix of another ("IntensityOccupancy_" does not start with "Intensity_"), and a
    // label is everything after its prefix, underscores included.
    private static readonly string[] SamplePrefixes =
        ["IntensityOccupancy_", "CountOccupancy_", "SpectralCount_", "Intensity_"];

    public override SupportedFileType FileType => SupportedFileType.MetaMorpheusQuantifiedProteinGroups;
    public override Software Software { get; set; }

    public ProteinGroupFromTsvFile(string filePath) : base(filePath, Software.MetaMorpheus) { }

    /// <summary>
    /// Constructor used to initialize from the factory method
    /// </summary>
    public ProteinGroupFromTsvFile() : base() { }

    /// <exception cref="MzLibException">The file is not a readable protein-group table.</exception>
    public override void LoadResults()
    {
        try
        {
            using var csv = new CsvReader(new StreamReader(FilePath), ProteinGroupFromTsv.CsvConfiguration);
            csv.Read();
            csv.ReadHeader();
            var sampleColumns = SampleColumns(csv.HeaderRecord ?? [], SamplePrefixes);

            var results = new List<ProteinGroupFromTsv>();
            while (csv.Read())
            {
                var row = csv.GetRecord<ProteinGroupFromTsv>();
                var groups = new Dictionary<string, SampleGroupMeasurement>();
                foreach (var (index, prefix, label) in sampleColumns)
                {
                    if (!groups.TryGetValue(label, out var m))
                        groups[label] = m = new SampleGroupMeasurement(label);
                    string cell = csv.GetField(index) ?? "";
                    switch (prefix)
                    {
                        case "SpectralCount_": m.SpectralCount = Blank(cell) ? null : int.Parse(cell, NumberStyles.Integer, CultureInfo.InvariantCulture); break;
                        case "Intensity_": m.Intensity = Blank(cell) ? null : double.Parse(cell, NumberStyles.Float, CultureInfo.InvariantCulture); break;
                        case "CountOccupancy_": m.CountOccupancyText = Blank(cell) ? null : cell; break;
                        case "IntensityOccupancy_": m.IntensityOccupancyText = Blank(cell) ? null : cell; break;
                    }
                }
                row.SampleGroups = groups;
                results.Add(row);
            }
            Results = results;
        }
        catch (Exception e) when (e is not MzLibException)
        {
            throw new MzLibException($"Could not read protein-group file '{FilePath}': {e.Message}", e);
        }
    }

    /// <exception cref="NotSupportedException">Always. The format is written by
    /// <c>Omics.BioPolymerGroup.BioPolymerGroupTsvSchema</c>, which holds the groups this file only describes.</exception>
    public override void WriteResults(string outputPath) =>
        throw new NotSupportedException(
            "Writing protein-group tables is not supported; they are written from BioPolymerGroups by BioPolymerGroupTsvSchema.");

    /// <summary>
    /// (column index, prefix, label) for every per-sample column, in header order. The first prefix a
    /// header starts with wins, so a prefix that another one extends must come after it. Shared with
    /// <see cref="QuantifiedPeptideFile"/>.
    /// </summary>
    internal static List<(int Index, string Prefix, string Label)> SampleColumns(string[] header, string[] prefixes)
    {
        var columns = new List<(int, string, string)>();
        for (int i = 0; i < header.Length; i++)
        {
            string? prefix = prefixes.FirstOrDefault(p => header[i].StartsWith(p, StringComparison.Ordinal));
            if (prefix != null)
                columns.Add((i, prefix, header[i][prefix.Length..]));
        }
        return columns;
    }

    private static bool Blank(string cell) => string.IsNullOrWhiteSpace(cell);
}
