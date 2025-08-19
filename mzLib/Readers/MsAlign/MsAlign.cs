using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System.Text.RegularExpressions;

namespace Readers;

/// <summary>
/// Parent class to define the shared behavior of MsAlign spectra file types
/// </summary>
public class MsAlign : MsDataFile
{

    #region Magic Numbers (defaults on parsing fail)

    private const int DefaultErrorValue = -1;
    private const double DefaultMinMz = 0;
    private const double DefaultMaxMz = 2000;
    private const double DefaultIsolationWidth = 3;
    private const int DefaultPrecursorCharge = 0;

    #endregion

    /// <summary>
    /// Ms2Align will sometimes have a header which states many parameters.
    /// Not all software output a parameters header in their MsAlign. 
    /// Not all of these are required for the MsAlign object, but they are stored here for reference.
    /// </summary>
    #region Optional MsAlign Header 

    private Dictionary<string, string> ParsedHeader { get; set; }
    
    public string? FileName { get; private set; }
    public bool? FaimsData { get; private set; }
    public double? FaimsVoltage { get; private set; }
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
    public bool? UseMsDeconvScore { get; private set; }

    #endregion

    /// <summary>
    /// An easy way to distinguish the different types of MsAlign files
    /// 0: Combined file (Exists only in our code)
    /// 1: ms1.msalign deconvolution output
    /// 2: ms2.msalign deconvolution output
    /// </summary>
    public virtual int DefaultMsnOrder { get; private set; }

    protected MsAlign(int numSpectra, SourceFile sourceFile) : base(numSpectra, sourceFile)
    {
        ParsedHeader = new();
    }

    protected MsAlign(MsDataScan[] scans, SourceFile sourceFile) : base(scans, sourceFile)
    {
        ParsedHeader = new();

        // Sort the scans by OneBasedScanNumber
        var orderedScans = scans.OrderBy(x => x.OneBasedScanNumber).ToArray();
        var indexedScans = new MsDataScan[orderedScans[^1].OneBasedScanNumber];
        foreach (var scan in orderedScans)
            indexedScans[scan.OneBasedScanNumber - 1] = scan;

        IndexedScans = indexedScans;
        Scans = orderedScans;
    }

    protected MsAlign(string filePath) : base(filePath)
    {
        ParsedHeader = new();
    }

    public override MsDataScan GetOneBasedScan(int scanNumber)
    {
        return IndexedScans[scanNumber - 1];
    }

    public override MsDataFile LoadAllStaticData(FilteringParams filteringParams = null, int maxThreads = 1)
    {
        List<MsDataScan> scans = new();

        // Instantiate enums to keep track of where we are in the file. The header is entirely optional. 
        var headerProgress = ReadingProgress.NotFound;
        var entryProgress = ReadingProgress.NotFound;
        using (var sr = new StreamReader(FilePath))
        {
            List<string> linesToProcess = new();
            while (sr.ReadLine() is { } line)
            {
                if (line.Contains("BEGIN IONS")) // Once we see this, we know we are past the header
                    headerProgress = ReadingProgress.Finished;

                // get header, line by line. 
                if (headerProgress != ReadingProgress.Finished)
                {
                    switch (headerProgress)
                    {
                        case ReadingProgress.NotFound when line.Contains("##### Parameters #####"): // we have found the beginning of the header
                            headerProgress = ReadingProgress.Found;
                            break;
                        case ReadingProgress.Found when line.Contains("##### Parameters #####"): //  we have reached the end of the header, process what we can into the ParsedHeader
                            headerProgress = ReadingProgress.Finished;
                            ParseHeaderLines(linesToProcess);
                            linesToProcess.Clear();
                            break;
                        default: // still in the middle of the header, keep collecting the lines
                            linesToProcess.Add(line);
                            continue;
                    }
                }
                else // read data scan entries
                {
                    switch (entryProgress)
                    {
                        // each entry after header
                        case ReadingProgress.NotFound when line.Contains("BEGIN IONS"): // we have found the beginning of a new entry
                            entryProgress = ReadingProgress.Found;
                            break;
                        case ReadingProgress.Found when line.Contains("END IONS"): // we have found the end of our current entry, process all lines into a NeutralMassSpectrum
                            {
                                entryProgress = ReadingProgress.NotFound;
                                var scan = ParseEntryLines(linesToProcess, filteringParams, PrecursorWindowSize);
                                scans.Add(scan);
                                linesToProcess.Clear();
                                break;
                            }
                        default: // we are in the middle of our current entry. 
                            linesToProcess.Add(line);
                            break;
                    }
                }
            }
        }

        SourceFile = GetSourceFile();

        // ensures that if a scan (OneBasedScanNumber) does not exist,
        // the final scans array will contain a null value  
        // this unique case is due to the nature of loading msalign files
        var orderedScans = scans.OrderBy(x => x.OneBasedScanNumber).ToArray();
        var indexedScans = new MsDataScan[orderedScans[^1].OneBasedScanNumber];
        foreach (var scan in orderedScans)
            indexedScans[scan.OneBasedScanNumber - 1] = scan;

        IndexedScans = indexedScans;
        Scans = orderedScans;
        Software ??= Readers.Software.Unspecified;
        return this;
    }

    private SourceFile? _sourceFile;

    public override SourceFile GetSourceFile()
    {
        if (_sourceFile is not null)
            return _sourceFile;

        string trimmedFilePath = FilePath.GetPeriodTolerantFilenameWithoutExtension()
            .Replace("_ms1", "")
            .Replace("_ms2", "");
        return _sourceFile = 
            new SourceFile("no nativeID format", $"ms{DefaultMsnOrder}.msalign format", null, null, trimmedFilePath, DefaultMsnOrder.ToString());
    }

    #region Dynamic Connection

    protected MsDataScan[] IndexedScans { get; set; }

    private StreamReader? _streamReader;
    private Dictionary<int, long> _scanByteOffset;
    private static readonly Regex ScanNumberParser = new Regex(@"(^|\s)SCANS=(.*?)($|\D)");

    public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null)
    {
        lock (DynamicReadingLock)
        {
            if (_streamReader == null)
            {
                throw new MzLibException("Cannot get scan; the dynamic data connection to " + FilePath + " has been closed!");
            }

            if (_scanByteOffset.TryGetValue(oneBasedScanNumber, out long byteOffset))
            {
                // seek to the byte of the scan
                _streamReader.BaseStream.Position = byteOffset;
                _streamReader.DiscardBufferedData();

                return GetNextMsDataOneBasedScanFromConnection(_streamReader, filterParams);
            }
            else
            {
                throw new MzLibException("The specified scan number: " + oneBasedScanNumber + " does not exist in " + FilePath);
            }
        }
    }

    public override void CloseDynamicConnection() => _streamReader?.Dispose();

    public override void InitiateDynamicConnection()
    {
        if (!File.Exists(FilePath))
            throw new FileNotFoundException();
        if (_streamReader is not null && _streamReader.BaseStream.CanRead)
            return;

        _streamReader = new StreamReader(FilePath);

        BuildIndex();
    }

    private void BuildIndex()
    {
        _scanByteOffset = new Dictionary<int, long>();
        int oneBasedScanNumber = 0;
        long oneBasedScanByteOffset = 0;
        bool scanHasAScanNumber = false;

        while (_streamReader != null && _streamReader.Peek() > 0)
        {
            var currentPositionByteOffset = TextFileReading.GetByteOffsetAtCurrentPosition(_streamReader);

            string? line = _streamReader.ReadLine();

            if (line != null && line.StartsWith("BEGIN IONS", StringComparison.InvariantCultureIgnoreCase))
            {
                oneBasedScanByteOffset = currentPositionByteOffset;
                scanHasAScanNumber = false;
            }
            else if (line != null && line.StartsWith("SCANS=", StringComparison.InvariantCultureIgnoreCase))
            {
                scanHasAScanNumber = true;

                Match result = ScanNumberParser.Match(line);
                var scanString = result.Groups[2].Value;
                oneBasedScanNumber = int.Parse(scanString);
            }
            else if (line != null && line.StartsWith("END IONS", StringComparison.InvariantCultureIgnoreCase))
            {
                if (!scanHasAScanNumber)
                {
                    oneBasedScanNumber++;
                }

                if (!_scanByteOffset.TryAdd(oneBasedScanNumber, oneBasedScanByteOffset))
                {
                    throw new MzLibException("Scan number " + oneBasedScanNumber.ToString() +
                                             " appeared multiple times in " + FilePath + ", which is not allowed because we assume all scan numbers are unique.");
                }
            }
        }
    }

    private MsDataScan GetNextMsDataOneBasedScanFromConnection(StreamReader sr, IFilteringParams filterParams = null)
    {
        var entryProgress = ReadingProgress.NotFound;
        List<string> linesToProcess = new();
        // read the scan data
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
                    goto FoundAllLines;
                }
                default:
                    linesToProcess.Add(line);
                    break;
            }
        }

        FoundAllLines:

        return ParseEntryLines(linesToProcess, filterParams, PrecursorWindowSize);
    }

    #endregion

    #region Scan Infomration Parsing

    private void ParseHeaderLines(List<string> headerLines)
    {
        ParsedHeader = new Dictionary<string, string>();

        foreach (var line in headerLines)
        {
            if (line.StartsWith("#"))
            {
                var keyValue = line.TrimStart('#').Split(':');
                if (keyValue.Length == 2)
                {
                    var key = keyValue[0].Trim();
                    var value = keyValue[1].Trim();
                    ParsedHeader[key] = value;
                }
                else if (line.Contains("IsoDec", StringComparison.InvariantCultureIgnoreCase))
                    Software = Readers.Software.IsoDec;
                else if (line.Contains("TopFD", StringComparison.InvariantCultureIgnoreCase))
                    Software = Readers.Software.TopFD;
            }
        }

        // Set the parsed header values to the corresponding properties
        foreach (var key in ParsedHeader.Keys)
        {
            switch (key)
            {
                case "File name":
                    FileName = ParsedHeader[key];
                    break;
                case "Faims data":
                    FaimsData = ParsedHeader[key] == "Yes";
                    break;
                case "Faims voltage":
                    FaimsVoltage = ParsedHeader[key]?.ToNullableDouble();
                    break;
                case "Activation type":
                    DissociationType = Enum.TryParse(ParsedHeader[key], out DissociationType dissType)
                        ? dissType
                        : MassSpectrometry.DissociationType.Autodetect;
                    break;
                case "Number of MS1 scans":
                    Ms1ScanCount = ParsedHeader[key]?.ToNullableInt();
                    break;
                case "Number of MS2 scans":
                case "Number of MS/MS scans":
                    Ms2ScanCount = ParsedHeader[key]?.ToNullableInt();
                    break;
                case "Spectral data type":
                    SpectralDataType = ParsedHeader[key];
                    break;
                case "Maximum charge":
                    MaxAssumedChargeState = ParsedHeader[key]?.ToNullableInt();
                    break;
                case "Maximum monoisotopic mass":
                    MaxAssumedMonoisotopicMass = ParsedHeader[key].Split(" ")[0].ToNullableDouble();
                    break;
                case "Peak error tolerance":
                    PeakErrorTolerance = ParsedHeader[key];
                    break;
                case "MS1 signal/noise ratio":
                    Ms1SnRRatio = ParsedHeader[key]?.ToNullableDouble();
                    break;
                case "MS/MS signal/noise ratio":
                    Ms2SnRRatio = ParsedHeader[key]?.ToNullableDouble();
                    break;
                case "Thread number":
                    MaxThreadsToUse = ParsedHeader[key]?.ToNullableDouble();
                    break;
                case "Default precursor window":
                case "Precursor window size":
                    PrecursorWindowSize = ParsedHeader[key].Replace("m/z", "").Trim().ToNullableDouble();
                    break;
                case "Use MS-Deconv score":
                    UseMsDeconvScore = ParsedHeader[key] == "Yes";
                    break;
                case "Use Env CNN model":
                    UseEnvCnnModel = ParsedHeader[key] == "Yes";
                    break;
                case "Miss MS1 spectra":
                    MissMs1Spectra = ParsedHeader[key] == "Yes";
                    break;
                case "Version":
                    SoftwareVersion = ParsedHeader[key];
                    break;
                default:
                    break;
            }
        }
    }

    private MsDataScan ParseEntryLines(List<string> entryLines, IFilteringParams? filteringParams = null,
        double? isolationWidth = 3)
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
        double? precursorMzStart = null;
        double? precursorMzEnd = null;

        // This switch has all scan header properties that I have seen, including those which are not currently added to the MsDataScan
        foreach (var headerLine in entryLines.Where(p => p.Contains('=')))
        {
            var splits = headerLine.Split('=');
            switch (splits[0])
            {
                case "ID":
                case "SPECTRUM ID":
                    id = int.TryParse(splits[1], out int idValue) ? idValue : DefaultErrorValue;
                    break;
                case "FRACTION_ID":
                    fractionId = int.TryParse(splits[1], out int fractionIdValue) ? fractionIdValue : DefaultErrorValue;
                    break;
                case "FILE_NAME":
                    fileName = splits[1];
                    break;
                case "SCANS":
                    oneBasedScanNumber = int.TryParse(splits[1], out int scanNumberValue) ? scanNumberValue : DefaultErrorValue;
                    break;
                case "RETENTION_TIME":
                    // Divide by sixty as msAligns are in seconds and MsDataScan is in minutes
                    retentionTime = double.TryParse(splits[1], out double retentionTimeValue) ? retentionTimeValue / 60.0 : DefaultErrorValue;
                    break;
                case "LEVEL":
                    msnOrder = int.TryParse(splits[1], out int msnOrderValue) ? msnOrderValue : DefaultMsnOrder;
                    break;
                case "ACTIVATION":
                    dissociationType = Enum.TryParse(splits[1], true, out DissociationType result)
                        ? result : MassSpectrometry.DissociationType.Autodetect;
                    break;
                case "MS_ONE_ID":
                    precursorScanId = int.TryParse(splits[1], out int scanIdValue) ? scanIdValue : DefaultErrorValue;
                    break;
                case "MS_ONE_SCAN":
                    oneBasedPrecursorScanNumber = int.TryParse(splits[1], out int precursorScanNumberValue) ? precursorScanNumberValue : DefaultErrorValue;
                    break;
                case "PRECURSOR_MZ":
                    precursorMz = double.TryParse(splits[1], out double mzValue) ? mzValue : DefaultErrorValue;
                    break;
                case "PRECURSOR_CHARGE":
                    precursorCharge = int.TryParse(splits[1], out int chargeValue) ? chargeValue : DefaultPrecursorCharge;
                    break;
                case "PRECURSOR_MASS":
                    precursorMass = double.TryParse(splits[1], out double massValue) ? massValue : DefaultErrorValue;
                    break;
                case "PRECURSOR_INTENSITY":
                    precursorIntensity = double.TryParse(splits[1], out double intensityValue) ? intensityValue : DefaultErrorValue;
                    break;
                case "PRECURSOR_WINDOW_BEGIN":
                    precursorMzStart = double.TryParse(splits[1], out double mzStartValue) ? mzStartValue : DefaultErrorValue;
                    break;
                case "PRECURSOR_WINDOW_END":
                    precursorMzEnd = double.TryParse(splits[1], out double mzEndValue) ? mzEndValue : DefaultErrorValue;
                    break;

            }
        }

        var peakLines = entryLines.Where(p => p.Contains('\t')).ToArray();
        var monoMasses = new double[peakLines.Length];
        var mzs = new double[peakLines.Length];
        var intensities = new double[peakLines.Length];
        var charges = new int[peakLines.Length];

        for (int i = 0; i < peakLines.Length; i++)
        {
            var splits = peakLines[i].Split('\t');

            charges[i] = int.Parse(splits[2]);
            monoMasses[i] = double.Parse(splits[0]);
            intensities[i] = double.Parse(splits[1]);
            mzs[i] = monoMasses[i].ToMz(charges[i]);
        }



        double minmz = mzs.Length == 0 ? DefaultMinMz : mzs.Min();
        double maxmz = mzs.Length == 0 ? DefaultMaxMz : mzs.Max();

        double? isolationMz = precursorMz;
        if (msnOrder == 1)
        {
            isolationWidth = null;
        }
        else if (precursorMzStart is not null && precursorMzEnd is not null)
        {
            isolationWidth = precursorMzEnd.Value - precursorMzStart.Value;
            isolationMz = precursorMzStart.Value + isolationWidth.Value / 2.0;
        }
        else
        {
            isolationWidth ??= DefaultIsolationWidth;
        }

        var sorted =
            monoMasses
                .Select((_, i) => new
                {
                    MonoMass = monoMasses[i],
                    Intensity = intensities[i],
                    Charge = charges[i]
                })
                .OrderBy(x => x.MonoMass)
                .ToArray();
        var spectrum = new NeutralMassSpectrum(
            sorted.Select(p => p.MonoMass).ToArray(),
            sorted.Select(p => p.Intensity).ToArray(),
            sorted.Select(p => p.Charge).ToArray(),
            true);

        var dataScan = new MsDataScan(spectrum, oneBasedScanNumber, msnOrder, true, Polarity.Positive, retentionTime,
            new MzRange(minmz - 1, maxmz + 1), null, MZAnalyzerType.Orbitrap,
            intensities.Sum(), null, null, $"scan={oneBasedScanNumber}",
            precursorMz, precursorCharge, precursorIntensity, isolationMz,
            isolationWidth, dissociationType, oneBasedPrecursorScanNumber, precursorMz);

        return dataScan;
    }

    #endregion

    #region Combine Ms1 and Ms2 Align

    // This region is for taking Ms1 and Ms2 align results of the same raw file and combining them to look like a standard LC/MS file with Ms2s following their respective Ms1
    // It is currently not used for any specific information, but might be useful for future work of database searches from external deconvolution results 

    public static bool TryCombineMsAlign(Ms1Align? ms1Align, Ms2Align? ms2Align, out MsAlign? combinedMsAlign)
    {
        if (ms1Align == null || ms2Align == null)
        {
            combinedMsAlign = null;
            return false;
        }
        

        combinedMsAlign = CombineMsAlign(ms1Align, ms2Align);
        return true;
        
    }

    internal static MsAlign CombineMsAlign(Ms1Align ms1Align, Ms2Align ms2Align)
    {
        var scans = new List<MsDataScan>();
        scans.AddRange(ms1Align.Scans);
        scans.AddRange(ms2Align.Scans);
        scans.Sort((scan1, scan2) => scan1.OneBasedScanNumber.CompareTo(scan2.OneBasedScanNumber)); // Sort the scans by scan number

        var sourceFile = ms1Align.GetSourceFile();
        var combinedFile =  new MsAlign(scans.ToArray(), sourceFile);

        //TODO: transfer over any important header properties
        
        combinedFile.DefaultMsnOrder = 0;

        return combinedFile;
    }

    #endregion

    /// <summary>
    /// Enum is required as there are several different ways an msAlign header information is written
    /// </summary>
    private enum ReadingProgress
    {
        NotFound,
        Found,
        Finished
    }
}