using MassSpectrometry;

namespace Readers;

public class SnipCreator
{
    public void Snipper(string origDataFile, int startScan, int endScan)
    {
        
        FilteringParams filter = new FilteringParams(200, 0.01, 1, null, false, false, true);
        var reader = MsDataFileReader.GetDataFile(origDataFile);
        reader.LoadAllStaticData(filter, 1);

        var scans = reader.GetAllScansList();
        var scansToKeep = scans.Where(x => x.OneBasedScanNumber >= startScan && x.OneBasedScanNumber <= endScan).ToList();

        List<(int oneBasedScanNumber, int? oneBasedPrecursorScanNumber)> scanNumbers = new List<(int oneBasedScanNumber, int? oneBasedPrecursorScanNumber)>();

        foreach (var scan in scansToKeep)
        {
            if (scan.OneBasedPrecursorScanNumber.HasValue && (scan.OneBasedPrecursorScanNumber.Value - startScan + 1) >= 0)
            {
                scanNumbers.Add((scan.OneBasedScanNumber, scan.OneBasedPrecursorScanNumber));
            }
        }

        Dictionary<int, int> scanNumberMap = new Dictionary<int, int>();

        foreach (var scanNumber in scanNumbers)
        {
            if (!scanNumberMap.ContainsKey(scanNumber.oneBasedScanNumber))
            {
                scanNumberMap.Add(scanNumber.oneBasedScanNumber, scanNumber.oneBasedScanNumber - startScan + 1);
            }
            if (scanNumber.oneBasedPrecursorScanNumber.HasValue && !scanNumberMap.ContainsKey(scanNumber.oneBasedPrecursorScanNumber.Value))
            {
                scanNumberMap.Add(scanNumber.oneBasedPrecursorScanNumber.Value, scanNumber.oneBasedPrecursorScanNumber.Value - startScan + 1);
            }
        }
        List<MsDataScan> scansForTheNewFile = new List<MsDataScan>();


        foreach (var scanNumber in scanNumbers)
        {
            MsDataScan scan = scansToKeep.First(x => x.OneBasedScanNumber == scanNumber.oneBasedScanNumber);

            MsDataScan newDataScan = new MsDataScan(
                scan.MassSpectrum,
                scanNumberMap[scan.OneBasedScanNumber],
                scan.MsnOrder,
                scan.IsCentroid,
                scan.Polarity,
                scan.RetentionTime,
                scan.ScanWindowRange,
                scan.ScanFilter,
                scan.MzAnalyzer,
                scan.TotalIonCurrent,
                scan.InjectionTime,
                scan.NoiseData,
                scan.NativeId.Replace(scan.OneBasedPrecursorScanNumber.ToString(), scanNumberMap[scan.OneBasedScanNumber].ToString()),
                scan.SelectedIonMZ,
                scan.SelectedIonChargeStateGuess,
                scan.SelectedIonIntensity,
                scan.IsolationMz,
                scan.IsolationWidth,
                scan.DissociationType,
                scanNumberMap[scan.OneBasedPrecursorScanNumber.Value],
                scan.SelectedIonMonoisotopicGuessMz,
                scan.HcdEnergy
            );
            scansForTheNewFile.Add(newDataScan);
        }

        string outPath = origDataFile.Replace(".mzML", "_snip.mzML");

        SourceFile sourceFile = new SourceFile(reader.SourceFile.NativeIdFormat,
            reader.SourceFile.MassSpectrometerFileFormat, reader.SourceFile.CheckSum, reader.SourceFile.Uri.ToString(),
            reader.SourceFile.FileName);

        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(scansForTheNewFile.ToArray(),sourceFile), outPath, false);
    }
}