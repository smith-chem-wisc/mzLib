using MassSpectrometry;

namespace Readers;

public class Ms2Align : MsAlign
{
        
    /// <summary>
    /// Ms2Align will sometimes have a header which states many parameters.
    /// Not all software output a parameters header in their Ms2Align. 
    /// Not all of these are required for the Ms2Align object, but they are stored here for reference.
    /// </summary>
    #region Optional Ms2Align Header 

    private Dictionary<string, string> _parsedHeader { get; set; }
    private void ParseHeaderLines(List<string> headerLines)
    {
        // TODO: Add header properties from IsoDec
        foreach (var line in headerLines.Where(p => p.Contains('\t')))
        {
            var splits = line.Split('\t');
            _parsedHeader.Add(splits[0].Trim(), splits[1].Trim());
        }

        // TODO: Parse the rest of the _parsedHeader Object into Header Properties
        string? precursorWindowSize = _parsedHeader.FirstOrDefault(p =>
            p.Key.Contains("precursor window", StringComparison.CurrentCultureIgnoreCase))!.Value;
        if (precursorWindowSize is not null) 
            PrecursorWindowSize = double.Parse(precursorWindowSize.Replace("m/z", "").Trim());
    }

    public DissociationType? DissociationType { get; private set; }
    public int? Ms1ScanCount { get; private set; }
    public int? Ms2ScanCount { get; private set; }
    public string? SpectralDataType { get; private set; }
    public int? MaxAssumedChargeState { get; private set; }
    public double? MaxAssumedMonoisotopicMass { get; private set; }
    public string? PeakErrorTolerance { get; private set; }
    public double? Ms1SnRRatio { get; private set; }
    public double? Ms2SnRRatio { get; private set; }
    public double? MaxThreadsToUse { get; private set; }
    public double? PrecursorWindowSize { get; private set; }
    public bool? UseEnvCnnModel { get; private set; }
    public bool? MissMs1Spectra { get; private set; }
    public string? SoftwareVersion { get; private set; }
    public Software? Software { get; private set; }

    #endregion


    protected override int DefaultMsnOrder => 2;

    public Ms2Align(int numSpectra, SourceFile sourceFile) : base(numSpectra, sourceFile)
    {
        _parsedHeader = [];
    }

    public Ms2Align(MsDataScan[] scans, SourceFile sourceFile) : base(scans, sourceFile)
    {
        _parsedHeader = [];
    }

    public Ms2Align(string filePath) : base(filePath) {
        _parsedHeader = [];
    }

    public override MsDataFile LoadAllStaticData(FilteringParams filteringParams = null, int maxThreads = 1)
    {
        List<MsDataScan> scans = [];
        var headerProgress = ReadingProgress.NotFound;
        var entryProgress = ReadingProgress.NotFound;
        using (var sr = new StreamReader(FilePath))
        {
            List<string> linesToProcess = [];
            while (sr.ReadLine() is { } line)
            {
                if (line.Contains("BEGIN IONS"))
                    headerProgress = ReadingProgress.Finished;

                // get header
                if (headerProgress != ReadingProgress.Finished)
                {
                    switch (headerProgress)
                    {
                        case ReadingProgress.NotFound when line.Contains("##### Parameters #####"):
                            headerProgress = ReadingProgress.Found;
                            break;
                        case ReadingProgress.Found when line.Contains("##### Parameters #####"):
                            headerProgress = ReadingProgress.Finished;
                            ParseHeaderLines(linesToProcess);
                            linesToProcess.Clear();
                            break;
                        default:
                            linesToProcess.Add(line);
                            continue;
                    }
                }
                else
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
                            var scan = PrecursorWindowSize is null
                                ? ParseEntryLines(linesToProcess, filteringParams)
                                : ParseEntryLines(linesToProcess, filteringParams, PrecursorWindowSize.Value);
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
        }

        Scans = scans.ToArray();
        return this;
    }

        
}