namespace Readers;

public class PsmFromTsvFile : SpectrumMatchFromTsvFile<PsmFromTsv>
{
    public override SupportedFileType FileType => SupportedFileType.psmtsv;
    public override Software Software { get; set; }

    /// <summary>
    /// Constructor used to initialize from the factory method
    /// </summary>
    public PsmFromTsvFile() : base() { }

    public PsmFromTsvFile(string filePath, SpectrumMatchParsingParameters? parsingParams = null) : base(filePath, parsingParams) { }

    public override void LoadResults()
    {
        Results = SpectrumMatchTsvReader.ReadPsmTsv(FilePath, out List<string> warnings, ParsingParams);
    }

    public override void WriteResults(string outputPath) => throw new NotImplementedException();
}