using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System.Text.RegularExpressions;

namespace Readers;

/// <summary>
/// Parent class to define the shared behavior of MsAlign spectra file types
/// </summary>
public abstract class MsAlign : MsDataFile
{
    // TODO: Add header properties from IsoDec and Ms1 Align from TopFD
    /// <summary>
    /// Ms2Align will sometimes have a header which states many parameters.
    /// Not all software output a parameters header in their MsAlign. 
    /// Not all of these are required for the MsAlign object, but they are stored here for reference.
    /// </summary>
    #region Optional MsAlign Header 

    protected Dictionary<string, string> _parsedHeader { get; set; }
   

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

    protected abstract int DefaultMsnOrder { get; }

    protected MsAlign(int numSpectra, SourceFile sourceFile) : base(numSpectra, sourceFile)
    {
        _parsedHeader = [];
    }

    protected MsAlign(MsDataScan[] scans, SourceFile sourceFile) : base(scans, sourceFile)
    {
        _parsedHeader = [];
    }

    protected MsAlign(string filePath) : base(filePath)
    {
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
    protected void ParseHeaderLines(List<string> headerLines)
    {
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
    protected MsDataScan ParseEntryLines(List<string> entryLines, IFilteringParams? filteringParams = null,
        double isolationWidth = 3)
    {
        // all
        int id;
        int fractionId;
        string fileName;
        int oneBasedScanNumber = 0;
        double retentionTime = 0;
        int msnOrder = 0;

        // ms2
        DissociationType? dissociationType = null;
        int? precursorScanId = null;
        int? oneBasedPrecursorScanNumber = null;
        double? precursorMz = null;
        int? precursorCharge = null;
        double? precursorMass = null;
        double? precursorIntensity = null;

        // This switch has all scan header properties that I have seen, including those which are not currently added to the MsDataScan
        foreach (var headerLine in entryLines.Where(p => p.Contains('=')))
        {
            var splits = headerLine.Split('=');
            switch (splits[0])
            {
                case "ID":
                case "SPECTRUM ID":
                    id = int.TryParse(splits[1], out int idValue) ? idValue : -1;
                    break;
                case "FRACTION_ID":
                    fractionId = int.TryParse(splits[1], out int fractionIdValue) ? fractionIdValue : -1;
                    break;
                case "FILE_NAME":
                    fileName = splits[1];
                    break;
                case "SCANS":
                    oneBasedScanNumber = int.TryParse(splits[1], out int scanNumberValue) ? scanNumberValue : -1;
                    break;
                case "RETENTION_TIME":
                    retentionTime = double.TryParse(splits[1], out double retentionTimeValue) ? retentionTimeValue : -1;
                    break;
                case "LEVEL":
                    msnOrder = int.TryParse(splits[1], out int msnOrderValue) ? msnOrderValue : DefaultMsnOrder;
                    break;
                case "ACTIVATION":
                    dissociationType = Enum.TryParse(splits[1], true, out DissociationType result)
                        ? result : MassSpectrometry.DissociationType.Autodetect;
                    break;
                case "MS_ONE_ID":
                    precursorScanId = int.TryParse(splits[1], out int scanIdValue) ? scanIdValue : -1;
                    break;
                case "MS_ONE_SCAN":
                    oneBasedPrecursorScanNumber = int.TryParse(splits[1], out int precursorScanNumberValue) ? precursorScanNumberValue : -1;
                    break;
                case "PRECURSOR_MZ":
                    precursorMz = double.TryParse(splits[1], out double mzValue) ? mzValue : -1;
                    break;
                case "PRECURSOR_CHARGE":
                    precursorCharge = int.TryParse(splits[1], out int chargeValue) ? chargeValue : 0;
                    break;
                case "PRECURSOR_MASS":
                    precursorMass = double.TryParse(splits[1], out double massValue) ? massValue : -1;
                    break;
                case "PRECURSOR_INTENSITY":
                    precursorIntensity = double.TryParse(splits[1], out double intensityValue) ? intensityValue : -1;
                    break;
            }
        }

        var peakLines = entryLines.Where(p => p.Contains('\t')).ToArray();
        var mzs = new double[peakLines.Length];
        var intensities = new double[peakLines.Length];
        var charges = new int[peakLines.Length];

        for (int i = 0; i < peakLines.Length; i++)
        {
            var splits = peakLines[i].Split('\t');

            charges[i] = int.Parse(splits[2]);
            mzs[i] = double.Parse(splits[0]).ToMz(charges[i]);
            intensities[i] = double.Parse(splits[1]);
        }

        double minMz = mzs.Min();
        double maxMz = mzs.Max();

        // peak filtering
        if (filteringParams != null && intensities.Length > 0 &&
            ((filteringParams.ApplyTrimmingToMs1 && msnOrder == 1) || (filteringParams.ApplyTrimmingToMsMs && msnOrder > 1)))
        {
            WindowModeHelper.Run(ref intensities, ref mzs, filteringParams, minMz, maxMz);
        }

        var spectrum = new MzSpectrum(mzs, intensities, true);
        var t = new MsDataScan(spectrum, oneBasedScanNumber, msnOrder, true, Polarity.Positive, retentionTime,
            mzs.Any() ? new MzRange(minMz, maxMz) : new MzRange(0, 2000), null, MZAnalyzerType.Orbitrap,
            intensities.Sum(), null, null, $"scan={oneBasedScanNumber}", precursorMz,
            precursorCharge, precursorIntensity, precursorMz, isolationWidth,
            dissociationType, oneBasedPrecursorScanNumber, precursorMass);

        //_parsedHeader.TryGetValue("Precursor window size:", out string value) ? double.Parse(value) : 3
        return t;
    }


    // TODO: Dynamic connection for MsAlign Types
    // Current approach is to have the dynamic methods call the static methods. 
    #region Dynamic Connection

    private StreamReader _streamReader;
    private Dictionary<int, long> _scanByteOffset;
    private static Regex _scanNumberparser = new Regex(@"(^|\s)SCANS=(.*?)($|\D)");

    public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null)
    {
        // TODO: Apply the filtering params
        return GetOneBasedScan(oneBasedScanNumber);
    }

    public override void CloseDynamicConnection() { }

    public override void InitiateDynamicConnection()
    {
        if (!CheckIfScansLoaded())
            LoadAllStaticData();
    }

    #endregion

    public override SourceFile GetSourceFile()
    {
        return new SourceFile("no nativeID format", $"ms{DefaultMsnOrder}.msalign format", null, null, null);
    }

    /// <summary>
    /// Enum is required as there are several different ways an msAlign header information is written
    /// </summary>
    protected enum ReadingProgress
    {
        NotFound,
        Found,
        Finished
    }
}