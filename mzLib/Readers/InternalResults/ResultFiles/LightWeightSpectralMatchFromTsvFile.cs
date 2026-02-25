namespace Readers;

/// <summary>
/// A result file wrapper for <see cref="LightWeightSpectralMatch"/> records.
/// Supports optional row and terminating filters that are passed through to
/// <see cref="LightWeightSpectralMatchTsvReader.ReadTsv"/>.
/// </summary>
public class LightWeightSpectralMatchFromTsvFile : ResultFile<LightWeightSpectralMatch>, IQuantifiableResultFile
{
    public override SupportedFileType FileType => SupportedFileType.psmtsv;
    public override Software Software { get; set; }

    /// <summary>
    /// Optional filters that skip non-matching rows.
    /// Key = column name from <see cref="SpectrumMatchFromTsvHeader"/>.
    /// Value = function taking the raw cell string, returning true to include the row.
    /// </summary>
    public Dictionary<string, Func<string, bool>>? RowFilters { get; set; }

    /// <summary>
    /// Optional filters that stop reading on first failure.
    /// Caller must ensure the file is sorted by the relevant column.
    /// </summary>
    public Dictionary<string, Func<string, bool>>? TerminatingFilters { get; set; }

    public LightWeightSpectralMatchFromTsvFile() : base() { }

    public LightWeightSpectralMatchFromTsvFile(string filePath) : base(filePath, Software.MetaMorpheus) { }

    public override void LoadResults()
    {
        Results = LightWeightSpectralMatchTsvReader.ReadTsv(FilePath, out _, RowFilters, TerminatingFilters);
    }

    public override void WriteResults(string outputPath) => throw new NotImplementedException();

    public IEnumerable<IQuantifiableRecord> GetQuantifiableResults() => Results;

    public Dictionary<string, string> FileNameToFilePath(List<string> fullFilePath)
    {
        var fileNames = Results.Select(r => r.FileName).Distinct().ToList();
        var fileNameToPath = new Dictionary<string, string>();
        foreach (string fullPath in fullFilePath)
        {
            var fileName = Path.GetFileNameWithoutExtension(fullPath);
            if (fileNames.Contains(fileName))
            {
                fileNameToPath[fileName] = fullPath;
            }
        }
        return fileNameToPath;
    }
}
