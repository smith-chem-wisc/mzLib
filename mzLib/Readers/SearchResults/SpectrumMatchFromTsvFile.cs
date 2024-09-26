using Omics.SpectrumMatch;

namespace Readers;

public class SpectrumMatchFromTsvFile : ResultFile<SpectrumMatchFromTsv>, IResultFile
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
}