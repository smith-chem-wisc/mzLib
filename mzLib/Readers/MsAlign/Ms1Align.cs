using MassSpectrometry;
namespace Readers;

public class Ms1Align : MsAlign
{
    protected override int DefaultMsnOrder => 1;

    public Ms1Align(int numSpectra, SourceFile sourceFile) : base(numSpectra, sourceFile) { }

    public Ms1Align(MsDataScan[] scans, SourceFile sourceFile) : base(scans, sourceFile) { }

    public Ms1Align(string filePath) : base(filePath) { }

    public override MsDataFile LoadAllStaticData(FilteringParams? filteringParams = null, int maxThreads = 1)
    {
        List<MsDataScan> scans = [];
        var entryProgress = ReadingProgress.NotFound;
        using (var sr = new StreamReader(FilePath))
        {
            List<string> linesToProcess = [];
            while (sr.ReadLine() is { } line)
            {
                switch (entryProgress)
                {
                    // each entry after header
                    case ReadingProgress.NotFound when line.Contains("BEGIN IONS"):
                        entryProgress = ReadingProgress.Found;
                        break;
                    case ReadingProgress.Found when line.Contains("END IONS"):
                    {
                        entryProgress = ReadingProgress.NotFound;
                        var scan = ParseEntryLines(linesToProcess, filteringParams);
                        scans.Add(scan);
                        linesToProcess.Clear();
                        break;
                    }
                    default:
                        linesToProcess.Add(line);
                        break;
                }
            }
        }

        Scans = scans.ToArray();
        return this;
    }
}