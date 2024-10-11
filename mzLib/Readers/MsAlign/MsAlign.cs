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

    public override bool IsNeutralMassFile { get; protected set; } = true;
    public abstract int DefaultMsnOrder { get; }

    protected MsAlign(int numSpectra, SourceFile sourceFile) : base(numSpectra, sourceFile)
    {
        ParsedHeader = new();
    }

    protected MsAlign(MsDataScan[] scans, SourceFile sourceFile) : base(scans, sourceFile)
    {
        ParsedHeader = new();
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
        var headerProgress = ReadingProgress.NotFound;
        var entryProgress = ReadingProgress.NotFound;
        using (var sr = new StreamReader(FilePath))
        {
            List<string> linesToProcess = new();
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
                                var scan = ParseEntryLines(linesToProcess, filteringParams, PrecursorWindowSize);
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

        SourceFile = GetSourceFile();

        // ensures that if a scan (OneBasedScanNumber) does not exist,
        // the final scans array will contain a null value  
        // this unique case is due to the nature of loading MGF files
        var orderedScans = scans.OrderBy(x => x.OneBasedScanNumber).ToArray();
        var indexedScans = new MsDataScan[orderedScans[^1].OneBasedScanNumber];
        foreach (var scan in orderedScans)
            indexedScans[scan.OneBasedScanNumber - 1] = scan;

        IndexedScans = indexedScans;
        Scans = orderedScans;
        Software ??= Readers.Software.Unspecified;
        return this;
    }

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
                    // Divide by sixty as msAligns are in seconds and MsDataScan is in minutes
                    retentionTime = double.TryParse(splits[1], out double retentionTimeValue) ? retentionTimeValue / 60.0 : -1;
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
                case "PRECURSOR_WINDOW_BEGIN":
                    precursorMzStart = double.TryParse(splits[1], out double mzStartValue) ? mzStartValue : -1;
                    break;
                case "PRECURSOR_WINDOW_END":
                    precursorMzEnd = double.TryParse(splits[1], out double mzEndValue) ? mzEndValue : -1;
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

        Array.Sort(mzs, intensities);

        //Remove Zero Intensity Peaks
        double zeroEquivalentIntensity = 0.01;
        int zeroIntensityCount = intensities.Count(i => i < zeroEquivalentIntensity);
        int intensityValueCount = intensities.Count();
        if (zeroIntensityCount > 0 && zeroIntensityCount < intensityValueCount)
        {
            Array.Sort(intensities, mzs);
            double[] nonZeroIntensities = new double[intensityValueCount - zeroIntensityCount];
            double[] nonZeroMzs = new double[intensityValueCount - zeroIntensityCount];
            intensities = intensities.SubArray(zeroIntensityCount, intensityValueCount - zeroIntensityCount);
            mzs = mzs.SubArray(zeroIntensityCount, intensityValueCount - zeroIntensityCount);
            Array.Sort(mzs, intensities);
        }

        double minMz = mzs.Length == 0 ? 0 : mzs.Min();
        double maxMz = mzs.Length == 0 ? 2000 : mzs.Max();

        // peak filtering
        if (filteringParams != null && intensities.Length > 0 &&
            ((filteringParams.ApplyTrimmingToMs1 && msnOrder == 1) || (filteringParams.ApplyTrimmingToMsMs && msnOrder > 1)))
        {
            WindowModeHelper.Run(ref intensities, ref mzs, filteringParams, minMz, maxMz);
        }

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
            isolationWidth ??= 3;
        }
        
        var spectrum = new MzSpectrum(mzs, intensities, true);
        var dataScan = new MsDataScan(spectrum, oneBasedScanNumber, msnOrder, true, Polarity.Positive, retentionTime,
            mzs.Any() ? new MzRange(minMz, maxMz) : new MzRange(0, 2000), null, MZAnalyzerType.Orbitrap,
            intensities.Sum(), null, null, $"scan={oneBasedScanNumber}",
            precursorMz, precursorCharge, precursorIntensity, isolationMz, 
            isolationWidth, dissociationType, oneBasedPrecursorScanNumber, precursorMz);

        return dataScan;
    }

    #region Dynamic Connection

    protected MsDataScan[] IndexedScans { get; set; }

    private StreamReader? _streamReader;
    private Dictionary<int, long> _scanByteOffset;
    private static readonly Regex ScanNumberParser = new Regex(@"(^|\s)SCANS=(.*?)($|\D)");

    public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null)
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

    private SourceFile? _sourceFile;
    public override SourceFile GetSourceFile() => _sourceFile ??= new SourceFile("no nativeID format", $"ms{DefaultMsnOrder}.msalign format", null, null, null);

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