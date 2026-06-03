namespace Readers;

public abstract class SpectrumMatchFromTsvFile<T> : ResultFile<T>, IQuantifiableResultFile where T: SpectrumMatchFromTsv
{
    protected SpectrumMatchParsingParameters? ParsingParams { get; }
    public override SupportedFileType FileType => FilePath.ParseFileType();
    public override Software Software { get; set; }

    /// <summary>
    /// Constructor used to initialize from the factory method
    /// </summary>
    protected SpectrumMatchFromTsvFile(SpectrumMatchParsingParameters? parsingParams = null) : base()
    {
        ParsingParams = parsingParams;
    }

    protected SpectrumMatchFromTsvFile(string filePath, SpectrumMatchParsingParameters? parsingParams = null) : base(filePath, Software.MetaMorpheus)
    {
        ParsingParams = parsingParams;
    }

    public override void LoadResults()
    {
        Results = SpectrumMatchTsvReader.ReadTsv<T>(FilePath, out List<string> warnings, ParsingParams);
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

public class SpectrumMatchFromTsvFile : SpectrumMatchFromTsvFile<SpectrumMatchFromTsv> 
{

    /// <summary>
    /// Constructor used to initialize from the factory method
    /// </summary>
    public SpectrumMatchFromTsvFile() : base() { }

    public SpectrumMatchFromTsvFile(string filePath) : base(filePath) { }
}