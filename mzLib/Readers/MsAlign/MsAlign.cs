using Chemistry;
using MassSpectrometry;
using MzLibUtil;

namespace Readers;

public abstract class MsAlign : MsDataFile
{
    /// <summary>
    /// Enum is required as there are several different ways an msAlign header information is written
    /// </summary>
    protected enum ReadingProgress
    {
        NotFound,
        Found,
        Finished
    }

    protected abstract int DefaultMsnOrder { get; }
    protected MsAlign(int numSpectra, SourceFile sourceFile) : base(numSpectra, sourceFile) { }
    protected MsAlign(MsDataScan[] scans, SourceFile sourceFile) : base(scans, sourceFile) { }
    protected MsAlign(string filePath) : base(filePath) { }

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
    public override SourceFile GetSourceFile()
    {
        throw new NotImplementedException();
    }
}