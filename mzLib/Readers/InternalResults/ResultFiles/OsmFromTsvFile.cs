using Transcriptomics;

namespace Readers;

public class OsmFromTsvFile : SpectrumMatchFromTsvFile<OsmFromTsv>
{
    public override SupportedFileType FileType => SupportedFileType.osmtsv;
    public override Software Software { get; set; }

    /// <summary>
    /// Constructor used to initialize from the factory method
    /// </summary>
    public OsmFromTsvFile() : base() { }

    public OsmFromTsvFile(string filePath) : base(filePath) { }

    public override void LoadResults()
    {
        Results = SpectrumMatchTsvReader.ReadOsmTsv(FilePath, out List<string> warnings);
    }

    public override void WriteResults(string outputPath) => throw new NotImplementedException();
}