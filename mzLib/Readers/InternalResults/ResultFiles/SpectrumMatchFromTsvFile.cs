using Omics.SpectrumMatch;

namespace Readers;

public class SpectrumMatchFromTsvFile : ResultFile<SpectrumMatchFromTsv>, IQuantifiableResultFile
{
    public override SupportedFileType FileType => FilePath.ParseFileType();
    public override Software Software { get; set; }

    /// <summary>
    /// Constructor used to initialize from the factory method
    /// </summary>
    public SpectrumMatchFromTsvFile() : base() { }

    public SpectrumMatchFromTsvFile(string filePath) : base(filePath, Software.MetaMorpheus) { }

    public override void LoadResults()
    {
        Results = SpectrumMatchTsvReader.ReadTsv(FilePath, out List<string> warnings);
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